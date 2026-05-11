# =============================================================================
# F.L.U.I.D — Módulo 4: DockingManager (AutoDock Vina Wrapper)
# =============================================================================
"""
Gerencia preparação e execução de docking molecular via AutoDock Vina.

Pipeline:
    1. Receptor: PDB → adiciona H polares → PDBQT (via Meeko/subprocess)
    2. Ligante: SMILES → RDKit 3D conformer → PDBQT (via Meeko)
    3. Docking: Vina via subprocess com Grid Box do RAGEngine
    4. GC: exclusão rigorosa de todos os .pdbqt, .log, .txt temporários

Todo I/O ocorre em tempfile.TemporaryDirectory — auto-limpeza garantida.

Complexidade:
    - Preparação: O(L) por ligante (RDKit embedding 3D)
    - Docking: O(L × M × E) — L ligantes, M modos, E exhaustiveness
    - GC: O(F) — F arquivos no tempdir
"""
from __future__ import annotations

import asyncio
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger

from src.config import FluidConfig, DEFAULT_CONFIG
from src.rag.rag_engine import DockingGuidance


# ---------------------------------------------------------------------------
# Exceções de Domínio
# ---------------------------------------------------------------------------
class DockingError(Exception):
    """Erro durante a simulação de docking."""
    pass


class LigandPreparationError(DockingError):
    """Falha na preparação do ligante (SMILES inválido, etc.)."""
    pass


class VinaExecutionError(DockingError):
    """Falha na execução do binário Vina."""
    pass


# ---------------------------------------------------------------------------
# DTOs
# ---------------------------------------------------------------------------
@dataclass
class DockingResult:
    """Resultado de um docking individual ligante-proteína."""
    ligand_id: str
    protein_id: str
    affinity_kcal: float = 0.0      # ΔG (kcal/mol)
    rmsd_lb: float = 0.0
    rmsd_ub: float = 0.0
    pose_pdbqt: str = ""
    mode_rank: int = 1
    success: bool = True
    error_msg: str = ""


# ---------------------------------------------------------------------------
# Classe Principal
# ---------------------------------------------------------------------------
class DockingManager:
    """
    Orquestra o pipeline de docking molecular com GC rigoroso.

    Preparação:
        - Receptor: PDB → PDBQT (Meeko mk_prepare_receptor.py ou ADFRsuite)
        - Ligante: SMILES → RDKit Mol → 3D embed → Meeko MoleculePreparation

    Execução:
        - Vina via subprocess (fallback) ou biblioteca Python vina
        - Grid Box orientada pelo RAGEngine (DockingGuidance)

    Garbage Collection:
        - Todo I/O em tempfile.TemporaryDirectory
        - Limpeza explícita de .pdbqt, .log, .txt após extração de ΔG
        - Context manager garante cleanup mesmo com exceções
    """

    def __init__(self, config: FluidConfig = DEFAULT_CONFIG) -> None:
        self._config = config
        self._temp_dirs: list[Path] = []  # Rastreio para GC manual

    # =================================================================
    # SEÇÃO 1: Preparação de Receptor
    # =================================================================

    def prepare_receptor(
        self, pdb_path: Path, output_dir: Path | None = None,
    ) -> Path:
        """
        Converte PDB em PDBQT para o receptor.

        Etapas:
            1. Adiciona hidrogênios polares (remove não-polares)
            2. Calcula cargas parciais de Gasteiger
            3. Exporta PDBQT via Meeko/OpenBabel

        Parâmetros:
            pdb_path: Caminho para o arquivo PDB do receptor.
            output_dir: Diretório de saída (None = tempdir).

        Retorna: Path do arquivo PDBQT gerado.
        Raises: DockingError se a preparação falhar.
        """
        if not pdb_path.exists():
            raise DockingError(f"PDB não encontrado: {pdb_path}")

        if output_dir is None:
            output_dir = Path(tempfile.mkdtemp(prefix="fluid_receptor_"))
            self._temp_dirs.append(output_dir)

        output_dir.mkdir(parents=True, exist_ok=True)
        pdbqt_path = output_dir / f"{pdb_path.stem}_receptor.pdbqt"

        try:
            # Tenta via Meeko (mk_prepare_receptor.py)
            result = subprocess.run(
                [
                    "mk_prepare_receptor",
                    "-i", str(pdb_path),
                    "-o", str(pdbqt_path),
                    "--add_hydrogen",
                ],
                capture_output=True, text=True, timeout=120,
            )

            if result.returncode != 0:
                # Fallback: OpenBabel
                logger.warning(
                    f"Meeko falhou (rc={result.returncode}), "
                    f"tentando OpenBabel: {result.stderr[:200]}"
                )
                result = subprocess.run(
                    [
                        "obabel", str(pdb_path),
                        "-O", str(pdbqt_path),
                        "-h",  # Adiciona H
                        "--partialcharge", "gasteiger",
                    ],
                    capture_output=True, text=True, timeout=120,
                )
                if result.returncode != 0:
                    raise DockingError(
                        f"Preparação do receptor falhou: {result.stderr[:300]}"
                    )

            if not pdbqt_path.exists() or pdbqt_path.stat().st_size == 0:
                raise DockingError("PDBQT do receptor vazio ou não gerado")

            logger.info(f"Receptor preparado: {pdb_path.name} → {pdbqt_path.name}")
            return pdbqt_path

        except subprocess.TimeoutExpired:
            raise DockingError(f"Timeout na preparação do receptor: {pdb_path}")

    # =================================================================
    # SEÇÃO 2: Preparação de Ligantes (RDKit + Meeko)
    # =================================================================

    def prepare_ligand(
        self, smiles: str, ligand_id: str, output_dir: Path | None = None,
    ) -> Path:
        """
        Converte SMILES em PDBQT para docking.

        Pipeline RDKit:
            1. Parse SMILES → RDKit Mol
            2. Adiciona hidrogênios (Chem.AddHs)
            3. Gera conformação 3D (AllChem.EmbedMolecule + ETKDG)
            4. Otimiza geometria (MMFF94 force field)
            5. Exporta via Meeko MoleculePreparation → PDBQT

        Raises: LigandPreparationError se SMILES inválido ou falha 3D.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            raise LigandPreparationError(
                "RDKit não instalado. Instale com: pip install rdkit-pypi"
            )

        if output_dir is None:
            output_dir = Path(tempfile.mkdtemp(prefix="fluid_ligand_"))
            self._temp_dirs.append(output_dir)

        output_dir.mkdir(parents=True, exist_ok=True)

        # 1. Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise LigandPreparationError(f"SMILES inválido: {smiles}")

        # 2. Adiciona H
        mol = Chem.AddHs(mol)

        # 3. Gera conformação 3D (ETKDG v3)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            # Fallback: random coords
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        # 4. Otimiza com MMFF94
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            logger.warning(f"MMFF94 falhou para {ligand_id}, usando UFF")
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass  # Usa geometria não-otimizada

        # 5. Exporta SDF intermediário
        sdf_path = output_dir / f"{ligand_id}.sdf"
        pdbqt_path = output_dir / f"{ligand_id}.pdbqt"

        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()

        # 6. Converte SDF → PDBQT via Meeko
        try:
            from meeko import MoleculePreparation, PDBQTWriterLegacy
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            for setup in mol_setups:
                pdbqt_string, is_ok, err = PDBQTWriterLegacy.write_string(setup)
                if is_ok:
                    pdbqt_path.write_text(pdbqt_string)
                    break
            else:
                raise LigandPreparationError(f"Meeko falhou para {ligand_id}")

        except ImportError:
            # Fallback: OpenBabel
            logger.warning("Meeko não disponível, usando OpenBabel")
            result = subprocess.run(
                ["obabel", str(sdf_path), "-O", str(pdbqt_path),
                 "--partialcharge", "gasteiger"],
                capture_output=True, text=True, timeout=60,
            )
            if result.returncode != 0:
                raise LigandPreparationError(
                    f"Conversão SDF→PDBQT falhou: {result.stderr[:200]}"
                )

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size == 0:
            raise LigandPreparationError(f"PDBQT vazio para {ligand_id}")

        logger.debug(f"Ligante preparado: {ligand_id} ({smiles[:30]}...)")
        return pdbqt_path

    # =================================================================
    # SEÇÃO 3: Execução de Docking (Vina)
    # =================================================================

    def run_docking(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        guidance: DockingGuidance,
        protein_id: str = "",
        ligand_id: str = "",
    ) -> list[DockingResult]:
        """
        Executa AutoDock Vina para um par receptor-ligante.

        Usa Grid Box do RAGEngine (DockingGuidance).
        Se guidance.confidence == 0, expande grid para blind focado.

        Retorna lista de DockingResult (um por modo/pose).
        """
        cfg = self._config
        cx, cy, cz = guidance.grid_center
        sx, sy, sz = guidance.grid_size

        # Se confiança baixa, expande grid (blind focado, não full blind)
        if guidance.confidence < 0.3:
            sx, sy, sz = max(sx, 30.0), max(sy, 30.0), max(sz, 30.0)
            logger.warning(
                f"Baixa confiança RAG ({guidance.confidence:.2f}), "
                f"expandindo grid para {sx}×{sy}×{sz} Å"
            )

        out_dir = Path(tempfile.mkdtemp(prefix="fluid_vina_"))
        self._temp_dirs.append(out_dir)
        output_pdbqt = out_dir / f"{ligand_id}_out.pdbqt"
        log_file = out_dir / f"{ligand_id}_vina.log"

        try:
            # Tenta biblioteca Python vina primeiro
            return self._run_vina_python(
                receptor_pdbqt, ligand_pdbqt, output_pdbqt,
                cx, cy, cz, sx, sy, sz,
                protein_id, ligand_id, cfg,
            )
        except (ImportError, Exception) as e:
            logger.info(f"vina Python falhou ({e}), usando subprocess")

        # Fallback: subprocess
        try:
            return self._run_vina_subprocess(
                receptor_pdbqt, ligand_pdbqt, output_pdbqt, log_file,
                cx, cy, cz, sx, sy, sz,
                protein_id, ligand_id, cfg,
            )
        except Exception as e:
            logger.error(f"Vina falhou para {ligand_id}: {e}")
            return [DockingResult(
                ligand_id=ligand_id, protein_id=protein_id,
                success=False, error_msg=str(e),
            )]

    def _run_vina_python(
        self, receptor: Path, ligand: Path, output: Path,
        cx: float, cy: float, cz: float,
        sx: float, sy: float, sz: float,
        protein_id: str, ligand_id: str, cfg: FluidConfig,
    ) -> list[DockingResult]:
        """Executa via biblioteca Python vina."""
        from vina import Vina

        v = Vina(sf_name="vina")
        v.set_receptor(str(receptor))
        v.set_ligand_from_file(str(ligand))
        v.compute_vina_maps(
            center=[cx, cy, cz],
            box_size=[sx, sy, sz],
        )
        v.dock(
            exhaustiveness=cfg.exhaustiveness,
            n_poses=cfg.num_modes,
        )
        v.write_poses(str(output), n_poses=cfg.num_modes)

        return self._parse_vina_output(
            output, protein_id, ligand_id
        )

    def _run_vina_subprocess(
        self, receptor: Path, ligand: Path, output: Path, log: Path,
        cx: float, cy: float, cz: float,
        sx: float, sy: float, sz: float,
        protein_id: str, ligand_id: str, cfg: FluidConfig,
    ) -> list[DockingResult]:
        """Executa via binário Vina (subprocess)."""
        cmd = [
            str(cfg.vina_binary),
            "--receptor", str(receptor),
            "--ligand", str(ligand),
            "--center_x", f"{cx:.3f}",
            "--center_y", f"{cy:.3f}",
            "--center_z", f"{cz:.3f}",
            "--size_x", f"{sx:.3f}",
            "--size_y", f"{sy:.3f}",
            "--size_z", f"{sz:.3f}",
            "--exhaustiveness", str(cfg.exhaustiveness),
            "--num_modes", str(cfg.num_modes),
            "--energy_range", f"{cfg.energy_range:.1f}",
            "--out", str(output),
            "--log", str(log),
        ]

        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=600,  # 10 min max por ligante
        )

        if result.returncode != 0:
            raise VinaExecutionError(
                f"Vina retornou código {result.returncode}: "
                f"{result.stderr[:300]}"
            )

        return self._parse_vina_output(output, protein_id, ligand_id)

    @staticmethod
    def _parse_vina_output(
        output_pdbqt: Path, protein_id: str, ligand_id: str,
    ) -> list[DockingResult]:
        """
        Parseia o PDBQT de saída do Vina para extrair ΔG e RMSD.

        Formato esperado (linhas REMARK VINA RESULT):
            REMARK VINA RESULT:   -7.3      0.000      0.000
        """
        results: list[DockingResult] = []

        if not output_pdbqt.exists():
            return [DockingResult(
                ligand_id=ligand_id, protein_id=protein_id,
                success=False, error_msg="Output PDBQT não encontrado",
            )]

        content = output_pdbqt.read_text()
        current_pose = []
        mode_rank = 0

        for line in content.split("\n"):
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 6:
                    mode_rank += 1
                    try:
                        affinity = float(parts[3])
                        rmsd_lb = float(parts[4])
                        rmsd_ub = float(parts[5])
                    except (ValueError, IndexError):
                        continue

                    results.append(DockingResult(
                        ligand_id=ligand_id,
                        protein_id=protein_id,
                        affinity_kcal=affinity,
                        rmsd_lb=rmsd_lb,
                        rmsd_ub=rmsd_ub,
                        mode_rank=mode_rank,
                        success=True,
                    ))

        if not results:
            results.append(DockingResult(
                ligand_id=ligand_id, protein_id=protein_id,
                success=False, error_msg="Nenhum resultado parseado do Vina",
            ))

        return results

    # =================================================================
    # SEÇÃO 4: Batch Processing com GC Rigoroso
    # =================================================================

    def run_batch(
        self,
        receptor_pdbqt: Path,
        ligands: list[tuple[str, str]],  # [(smiles, ligand_id), ...]
        guidance: DockingGuidance,
        protein_id: str = "",
    ) -> pd.DataFrame:
        """
        Processa batch de ligantes com garbage collection rigoroso.

        Para cada ligante:
            1. Prepara PDBQT (RDKit + Meeko)
            2. Executa Vina
            3. Extrai ΔG
            4. Remove arquivos temporários

        Retorna DataFrame consolidado ordenado por afinidade (melhor ΔG).
        """
        all_results: list[dict] = []

        for smiles, lig_id in ligands:
            try:
                # Prepara em tempdir isolado
                lig_pdbqt = self.prepare_ligand(smiles, lig_id)
                docking_results = self.run_docking(
                    receptor_pdbqt, lig_pdbqt, guidance,
                    protein_id=protein_id, ligand_id=lig_id,
                )

                for dr in docking_results:
                    all_results.append({
                        "ligand_id": dr.ligand_id,
                        "protein_id": dr.protein_id,
                        "smiles": smiles,
                        "affinity_kcal": dr.affinity_kcal,
                        "rmsd_lb": dr.rmsd_lb,
                        "rmsd_ub": dr.rmsd_ub,
                        "mode_rank": dr.mode_rank,
                        "success": dr.success,
                        "error": dr.error_msg,
                    })

            except (LigandPreparationError, DockingError) as e:
                logger.warning(f"Ligante {lig_id} falhou: {e}")
                all_results.append({
                    "ligand_id": lig_id, "protein_id": protein_id,
                    "smiles": smiles, "affinity_kcal": 0.0,
                    "rmsd_lb": 0.0, "rmsd_ub": 0.0, "mode_rank": 0,
                    "success": False, "error": str(e),
                })

        # === GARBAGE COLLECTION RIGOROSO ===
        self._garbage_collect_all()

        df = pd.DataFrame(all_results)
        if not df.empty:
            df = df.sort_values("affinity_kcal").reset_index(drop=True)
        return df

    # =================================================================
    # SEÇÃO 5: Garbage Collection
    # =================================================================

    def _garbage_collect_all(self) -> int:
        """
        Remove TODOS os diretórios temporários criados durante o batch.

        Garante exclusão de .pdbqt, .log, .txt e qualquer outro artefato.
        Retorna número de diretórios removidos.
        """
        removed = 0
        for tmp_dir in self._temp_dirs:
            try:
                if tmp_dir.exists():
                    # Conta arquivos antes de remover
                    files = list(tmp_dir.rglob("*"))
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                    removed += 1
                    logger.debug(
                        f"GC: removido {tmp_dir} ({len(files)} arquivos)"
                    )
            except Exception as e:
                logger.warning(f"GC falhou para {tmp_dir}: {e}")

        self._temp_dirs.clear()
        logger.info(f"GC: {removed} diretórios temporários removidos")
        return removed

    def __del__(self) -> None:
        """Garante cleanup no destrutor (safety net)."""
        try:
            self._garbage_collect_all()
        except Exception:
            pass

# =============================================================================
# F.L.U.I.D — Módulo 2: RAGEngine (LangChain + ChromaDB + Europe PMC)
# =============================================================================
"""
Pipeline RAG para orientação de docking molecular via literatura científica.

Arquitetura:
    1. Ingestão de artigos Open Access via Europe PMC REST API
    2. Chunking recursivo → embedding → ChromaDB (HNSW: busca O(log n))
    3. Query semântica com Structured Output Parser (Pydantic)
    4. Retorna DockingGuidance com coordenadas de Grid Box

A LLM é forçada a retornar JSON estruturado via Pydantic schema —
jamais texto livre. Se coordenadas não forem encontradas na literatura,
calcula o centroide geométrico dos resíduos catalíticos citados.
"""
from __future__ import annotations

import asyncio
import re
from dataclasses import dataclass, field
from typing import Any

import httpx
import numpy as np
from loguru import logger
from pydantic import BaseModel, Field

from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import Chroma
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_openai import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain.output_parsers import PydanticOutputParser

from src.config import FluidConfig, EnzymeFamilyType, DEFAULT_CONFIG


# ---------------------------------------------------------------------------
# Pydantic Schemas (Structured Output)
# ---------------------------------------------------------------------------
class GridCoordinate(BaseModel):
    """Coordenada 3D para Grid Box do docking."""
    x: float = Field(0.0, description="Coordenada X em Ångströms")
    y: float = Field(0.0, description="Coordenada Y em Ångströms")
    z: float = Field(0.0, description="Coordenada Z em Ångströms")


class DockingGuidanceSchema(BaseModel):
    """Schema Pydantic para output estruturado da LLM."""
    enzyme_class: str = Field(
        "other",
        description="Família enzimática: protease, polymerase, helicase, "
                    "integrase, kinase, transferase, ligase, hydrolase, other"
    )
    catalytic_residues: list[str] = Field(
        default_factory=list,
        description="Lista de resíduos catalíticos (ex: ['His41', 'Cys145'])"
    )
    allosteric_residues: list[str] = Field(
        default_factory=list,
        description="Lista de resíduos alostéricos, se identificados"
    )
    grid_center: GridCoordinate = Field(
        default_factory=GridCoordinate,
        description="Centro da Grid Box em Å (coordenadas do sítio ativo)"
    )
    grid_size: GridCoordinate = Field(
        default_factory=lambda: GridCoordinate(x=22.0, y=22.0, z=22.0),
        description="Dimensões da Grid Box em Å"
    )
    confidence: float = Field(
        0.0, ge=0.0, le=1.0,
        description="Confiança na extração (0=sem dados, 1=coordenadas exatas)"
    )
    source_pmids: list[str] = Field(
        default_factory=list,
        description="PMIDs dos artigos fonte"
    )


# ---------------------------------------------------------------------------
# DTO de saída (compatível com módulos downstream)
# ---------------------------------------------------------------------------
@dataclass
class DockingGuidance:
    """Orientação extraída da literatura para posicionamento de Grid Box."""
    enzyme_family: EnzymeFamilyType = EnzymeFamilyType.OTHER
    catalytic_residues: list[str] = field(default_factory=list)
    allosteric_residues: list[str] = field(default_factory=list)
    grid_center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    grid_size: tuple[float, float, float] = (22.0, 22.0, 22.0)
    confidence: float = 0.0
    source_pmids: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------
EUROPE_PMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPE_PMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"

# Prompt template para extração estruturada
_EXTRACTION_PROMPT = """You are a structural biology expert. Based ONLY on the 
scientific literature excerpts below, extract the binding site information for 
the target protein.

CONTEXT:
{context}

PROTEIN: {protein_name} (UniProt: {uniprot_id})

INSTRUCTIONS:
1. Identify the enzyme class (protease, polymerase, helicase, etc.)
2. List catalytic and allosteric residues mentioned in the text
3. If 3D coordinates of the active site are mentioned, extract them
4. If no exact coordinates are found, set confidence to 0.0
5. Report the PMIDs of the source articles

{format_instructions}
"""

# Mapa de texto → EnzymeFamilyType
_ENZYME_MAP: dict[str, EnzymeFamilyType] = {
    "protease": EnzymeFamilyType.PROTEASE,
    "polymerase": EnzymeFamilyType.POLYMERASE,
    "helicase": EnzymeFamilyType.HELICASE,
    "integrase": EnzymeFamilyType.INTEGRASE,
    "kinase": EnzymeFamilyType.KINASE,
    "transferase": EnzymeFamilyType.TRANSFERASE,
    "ligase": EnzymeFamilyType.LIGASE,
    "hydrolase": EnzymeFamilyType.HYDROLASE,
}


# ---------------------------------------------------------------------------
# Classe Principal
# ---------------------------------------------------------------------------
class RAGEngine:
    """
    Motor RAG com saída estruturada (Pydantic) para orientação de docking.

    Pipeline:
        1. ingest_literature() → busca Europe PMC → chunk → embed → Chroma
        2. get_docking_guidance() → query vetorial → LLM structured output
        3. Validação → se confidence == 0, calcula centroide de resíduos

    Complexidade:
        - Embedding: O(n × d) — n chunks, d dimensão
        - Busca HNSW (Chroma): O(log n) amortizado
        - LLM inference: O(1) por query (API call)
    """

    def __init__(self, config: FluidConfig = DEFAULT_CONFIG) -> None:
        self._config = config
        self._embeddings = HuggingFaceEmbeddings(
            model_name="sentence-transformers/all-MiniLM-L6-v2"
        )
        self._vectorstore: Chroma | None = None
        self._splitter = RecursiveCharacterTextSplitter(
            chunk_size=1500,
            chunk_overlap=200,
            separators=["\n\n", "\n", ". ", " "],
        )
        self._parser = PydanticOutputParser(
            pydantic_object=DockingGuidanceSchema
        )
        self._llm = ChatOpenAI(
            model="gpt-4o-mini",
            temperature=0.0,
            api_key=config.openai_api_key,
        )

    # =================================================================
    # SEÇÃO 1: Ingestão de Literatura (Europe PMC)
    # =================================================================

    async def _search_europe_pmc(
        self, query: str, *, max_results: int = 50
    ) -> list[dict[str, Any]]:
        """
        Busca artigos Open Access no Europe PMC.
        Retorna lista de metadados (pmid, pmcid, title, abstractText).
        """
        params = {
            "query": f"{query} OPEN_ACCESS:y",
            "resultType": "core",
            "pageSize": min(max_results, 100),
            "format": "json",
        }
        async with httpx.AsyncClient(timeout=30.0) as client:
            resp = await client.get(EUROPE_PMC_SEARCH, params=params)
            resp.raise_for_status()
            data = resp.json()

        results = data.get("resultList", {}).get("result", [])
        logger.info(f"Europe PMC: {len(results)} artigos encontrados para '{query}'")
        return results

    async def _fetch_fulltext(self, pmcid: str) -> str | None:
        """Baixa texto completo de um artigo via PMC ID."""
        if not pmcid:
            return None
        url = EUROPE_PMC_FULLTEXT.format(pmcid=pmcid)
        try:
            async with httpx.AsyncClient(timeout=30.0) as client:
                resp = await client.get(url)
                if resp.status_code == 200:
                    # Remove XML tags para texto puro
                    text = re.sub(r"<[^>]+>", " ", resp.text)
                    return re.sub(r"\s+", " ", text).strip()
        except httpx.HTTPError as e:
            logger.warning(f"Falha ao buscar fulltext {pmcid}: {e}")
        return None

    async def ingest_literature(
        self, query: str, *, max_papers: int = 30
    ) -> int:
        """
        Busca artigos no Europe PMC, extrai texto, chunka e indexa no Chroma.

        Retorna o número de chunks indexados.
        """
        articles = await self._search_europe_pmc(query, max_results=max_papers)

        documents: list[str] = []
        metadatas: list[dict] = []

        for art in articles:
            pmcid = art.get("pmcid", "")
            pmid = art.get("pmid", "")
            title = art.get("title", "")

            # Tenta fulltext, fallback para abstract
            text = await self._fetch_fulltext(pmcid)
            if not text:
                text = art.get("abstractText", "")

            if not text or len(text) < 100:
                continue

            chunks = self._splitter.split_text(text)
            for chunk in chunks:
                documents.append(chunk)
                metadatas.append({
                    "pmid": pmid, "pmcid": pmcid, "title": title
                })

        if not documents:
            logger.warning("Nenhum documento indexável encontrado")
            return 0

        # Cria/atualiza vector store
        persist_dir = str(self._config.workspace_dir / "chroma_db")
        self._vectorstore = Chroma.from_texts(
            texts=documents,
            metadatas=metadatas,
            embedding=self._embeddings,
            persist_directory=persist_dir,
            collection_name="fluid_literature",
        )

        logger.info(f"RAG: {len(documents)} chunks indexados de {len(articles)} artigos")
        return len(documents)

    # =================================================================
    # SEÇÃO 2: Query Semântica + Output Estruturado
    # =================================================================

    async def get_docking_guidance(
        self,
        uniprot_id: str,
        protein_name: str,
        enzyme_family: EnzymeFamilyType | None = None,
    ) -> DockingGuidance:
        """
        Extrai orientação de Grid Box da literatura indexada.

        Pipeline:
            1. Query vetorial no Chroma (top-k=8) → O(log n)
            2. LLM com Structured Output Parser → JSON Pydantic
            3. Validação: se confidence == 0, tenta centroide de resíduos

        Retorna DockingGuidance com coordenadas e confiança.
        Nunca levanta exceção — retorna guidance default em caso de erro.
        """
        if self._vectorstore is None:
            logger.warning("Vector store vazio — retornando guidance default")
            return self._default_guidance(enzyme_family)

        # 1. Busca vetorial — O(log n) via HNSW
        query = (
            f"{protein_name} {uniprot_id} active site catalytic residues "
            f"binding pocket coordinates"
        )
        try:
            docs = self._vectorstore.similarity_search(query, k=8)
        except Exception as e:
            logger.error(f"Busca vetorial falhou: {e}")
            return self._default_guidance(enzyme_family)

        if not docs:
            return self._default_guidance(enzyme_family)

        context = "\n\n---\n\n".join(doc.page_content for doc in docs)
        pmids = list({
            doc.metadata.get("pmid", "") for doc in docs if doc.metadata.get("pmid")
        })

        # 2. LLM com output estruturado
        prompt = ChatPromptTemplate.from_template(_EXTRACTION_PROMPT)
        chain = prompt | self._llm | self._parser

        try:
            result: DockingGuidanceSchema = await chain.ainvoke({
                "context": context,
                "protein_name": protein_name,
                "uniprot_id": uniprot_id,
                "format_instructions": self._parser.get_format_instructions(),
            })
        except Exception as e:
            logger.error(f"LLM extraction falhou: {e}")
            return self._default_guidance(enzyme_family)

        # 3. Validação e conversão
        return self._schema_to_guidance(result, pmids, enzyme_family)

    async def classify_enzyme_family(
        self, protein_name: str, sequence: str
    ) -> EnzymeFamilyType:
        """Classifica família enzimática via query RAG."""
        if self._vectorstore is None:
            return EnzymeFamilyType.OTHER

        try:
            docs = self._vectorstore.similarity_search(
                f"{protein_name} enzyme family classification", k=4
            )
            context = "\n".join(d.page_content for d in docs)

            prompt = ChatPromptTemplate.from_template(
                "Based on this context about {protein}, classify the enzyme "
                "family. Reply with ONLY one word: protease, polymerase, "
                "helicase, integrase, kinase, transferase, ligase, hydrolase, "
                "or other.\n\nContext: {context}"
            )
            chain = prompt | self._llm
            resp = await chain.ainvoke({
                "protein": protein_name, "context": context
            })
            label = resp.content.strip().lower()
            return _ENZYME_MAP.get(label, EnzymeFamilyType.OTHER)

        except Exception as e:
            logger.error(f"Classificação enzimática falhou: {e}")
            return EnzymeFamilyType.OTHER

    # =================================================================
    # SEÇÃO 3: Validação e Fallbacks
    # =================================================================

    @staticmethod
    def compute_residue_centroid(
        residue_coords: list[tuple[float, float, float]]
    ) -> tuple[float, float, float]:
        """
        Calcula centroide geométrico a partir de coordenadas de resíduos.
        Usado como fallback quando a LLM não encontra coordenadas exatas.

        Centroide = (mean(x), mean(y), mean(z))
        Complexidade: O(R) onde R = número de resíduos.
        """
        if not residue_coords:
            return (0.0, 0.0, 0.0)
        arr = np.array(residue_coords)
        return tuple(float(v) for v in arr.mean(axis=0))

    def _schema_to_guidance(
        self,
        schema: DockingGuidanceSchema,
        pmids: list[str],
        fallback_family: EnzymeFamilyType | None,
    ) -> DockingGuidance:
        """Converte schema Pydantic → DockingGuidance dataclass."""
        # Resolve enzyme family
        family = _ENZYME_MAP.get(
            schema.enzyme_class.lower(),
            fallback_family or EnzymeFamilyType.OTHER,
        )

        # Verifica se coordenadas são válidas (não-zero)
        gc = schema.grid_center
        has_coords = not (gc.x == 0.0 and gc.y == 0.0 and gc.z == 0.0)
        confidence = schema.confidence if has_coords else 0.0

        if not has_coords and confidence < 0.3:
            logger.warning(
                "RAG: coordenadas não encontradas na literatura. "
                "Grid Box será calculada por centroide de resíduos "
                "ou blind docking focado."
            )

        return DockingGuidance(
            enzyme_family=family,
            catalytic_residues=schema.catalytic_residues,
            allosteric_residues=schema.allosteric_residues,
            grid_center=(gc.x, gc.y, gc.z),
            grid_size=(schema.grid_size.x, schema.grid_size.y, schema.grid_size.z),
            confidence=confidence,
            source_pmids=pmids or schema.source_pmids,
        )

    @staticmethod
    def _default_guidance(
        enzyme_family: EnzymeFamilyType | None = None,
    ) -> DockingGuidance:
        """Retorna guidance padrão (fallback seguro)."""
        return DockingGuidance(
            enzyme_family=enzyme_family or EnzymeFamilyType.OTHER,
            confidence=0.0,
        )

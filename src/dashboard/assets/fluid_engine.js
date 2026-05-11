/**
 * =============================================================================
 * F.L.U.I.D — Motor de Física 2D: Simulação de Docking em Fluido
 * Engine: matter.js (física + broadphase grid) + p5.js (renderização)
 * =============================================================================
 *
 * Arquitetura de Complexidade:
 *   - Detecção de colisão: O(n) amortizado via Broadphase Grid (Spatial Hashing)
 *     O matter.js usa internamente Grid-based broadphase que particiona o espaço
 *     em células. Colisões são verificadas APENAS entre corpos na mesma célula
 *     + vizinhas (8-connectivity), reduzindo O(n²) → O(n).
 *
 *   - Movimento browniano: O(n) — uma força aleatória por partícula por frame
 *   - Atração por afinidade: O(m) — m = matches com ΔG < threshold
 *   - Renderização p5.js: O(n) — um draw call por corpo
 *
 * Fases da Simulação:
 *   1. FLUTUAÇÃO — Browniano com amortecimento viscoso (citosol)
 *   2. ATRAÇÃO — Constraint (mola) invisível ligante→proteína (ΔG < -7.0)
 *   3. ENCAIXE — Lock visual com animação de escala + opacidade
 */

// =============================================================================
// CONFIGURAÇÃO GLOBAL
// =============================================================================
const CONFIG = {
    // Canvas
    WIDTH: 900,
    HEIGHT: 600,
    BG_COLOR: '#0a0e1a',

    // Física
    BROWNIAN_FORCE: 0.0008,
    VISCOSITY: 0.05,           // frictionAir do matter.js
    AFFINITY_THRESHOLD: -7.0,  // kcal/mol — ΔG para ativar atração
    SPRING_STIFFNESS: 0.0004,  // Rigidez da mola de atração
    LOCK_DISTANCE: 15,         // px — distância para "encaixe"

    // Broadphase Grid — tamanho de célula para spatial hashing
    // Deve ser ≥ maior dimensão de corpo para eficiência O(n)
    GRID_BUCKET_SIZE: 80,

    // Visual
    PROTEIN_COLOR: '#4fc3f7',
    PROTEIN_GLOW: '#0288d1',
    LIGAND_COLORS: ['#ff7043', '#ab47bc', '#66bb6a', '#ffa726', '#ef5350'],
    MATCH_COLOR: '#76ff03',
    FONT_SIZE: 10,
};

// =============================================================================
// ESTADO GLOBAL
// =============================================================================
let engine, world, render_p5;
let proteins = [];
let ligands = [];
let constraints = [];
let matchedPairs = new Set();  // "ligandId::proteinId"
let dockingResults = {};       // { ligandId: { affinity, proteinId } }

// Matter.js aliases
const Engine = Matter.Engine;
const World = Matter.World;
const Bodies = Matter.Bodies;
const Body = Matter.Body;
const Constraint = Matter.Constraint;
const Events = Matter.Events;
const Vector = Matter.Vector;

// =============================================================================
// INICIALIZAÇÃO DO MOTOR DE FÍSICA
// =============================================================================

/**
 * Inicializa o motor matter.js com Broadphase Grid configurado.
 *
 * O broadphase.bucketWidth e bucketHeight definem o tamanho da célula
 * do spatial hash. Valores menores = mais células = menos comparações
 * por célula, mas mais overhead de hash. O sweet spot é ≈ 1.5× o
 * diâmetro médio dos corpos.
 */
function initPhysics() {
    engine = Engine.create({
        gravity: { x: 0, y: 0 },  // Sem gravidade (fluido 2D)
    });

    // Configura Broadphase Grid explicitamente
    // O matter.js ≥0.19 usa Grid por padrão, mas forçamos os parâmetros
    if (engine.broadphase) {
        engine.broadphase.bucketWidth = CONFIG.GRID_BUCKET_SIZE;
        engine.broadphase.bucketHeight = CONFIG.GRID_BUCKET_SIZE;
    }

    world = engine.world;

    // Paredes invisíveis (contêm as partículas)
    const wallOpts = { isStatic: true, render: { visible: false } };
    const w = CONFIG.WIDTH, h = CONFIG.HEIGHT, t = 30;
    World.add(world, [
        Bodies.rectangle(w / 2, -t / 2, w + 2 * t, t, wallOpts),
        Bodies.rectangle(w / 2, h + t / 2, w + 2 * t, t, wallOpts),
        Bodies.rectangle(-t / 2, h / 2, t, h + 2 * t, wallOpts),
        Bodies.rectangle(w + t / 2, h / 2, t, h + 2 * t, wallOpts),
    ]);

    // Listener de colisão para log (debug)
    Events.on(engine, 'collisionStart', (event) => {
        // Colisões gerenciadas pelo broadphase grid — O(n)
    });
}

// =============================================================================
// CRIAÇÃO DE CORPOS
// =============================================================================

/**
 * Cria um corpo de proteína (receptor) — grande, pesado, semi-estático.
 *
 * Representada como polígono convexo com "fenda" visual para o sítio ativo.
 * A forma é um hexágono irregular para simular a superfície proteica.
 */
function createProtein(id, label, x, y, size) {
    const vertices = [];
    const sides = 8;
    for (let i = 0; i < sides; i++) {
        const angle = (i / sides) * Math.PI * 2;
        // Irregularidade controlada para simular superfície proteica
        const r = size * (0.8 + 0.2 * Math.sin(i * 2.7));
        vertices.push({ x: r * Math.cos(angle), y: r * Math.sin(angle) });
    }

    const body = Bodies.fromVertices(x, y, vertices, {
        isStatic: false,
        mass: 50,         // Pesada — resiste ao browniano
        frictionAir: 0.1, // Alta viscosidade (quase estática)
        restitution: 0.3,
        label: `protein_${id}`,
    });

    body.fluidData = {
        type: 'protein',
        id: id,
        displayLabel: label,
        size: size,
        color: CONFIG.PROTEIN_COLOR,
        glowIntensity: 0,
        hasDocked: false,
    };

    World.add(world, body);
    proteins.push(body);
    return body;
}

/**
 * Cria um corpo de ligante — pequeno, leve, sujeito a browniano.
 *
 * Representado como círculo com sprite SVG do RDKit (quando disponível).
 */
function createLigand(id, label, svgData, x, y, radius) {
    radius = radius || 12 + Math.random() * 8;
    x = x || 50 + Math.random() * (CONFIG.WIDTH - 100);
    y = y || 50 + Math.random() * (CONFIG.HEIGHT - 100);

    const body = Bodies.circle(x, y, radius, {
        mass: 0.5,
        frictionAir: CONFIG.VISCOSITY,
        restitution: 0.6,
        label: `ligand_${id}`,
    });

    body.fluidData = {
        type: 'ligand',
        id: id,
        displayLabel: label,
        svgData: svgData || null,
        radius: radius,
        color: CONFIG.LIGAND_COLORS[ligands.length % CONFIG.LIGAND_COLORS.length],
        phase: 'floating',  // 'floating' | 'attracting' | 'docked'
        targetProtein: null,
        affinity: 0,
        opacity: 1.0,
        scale: 1.0,
    };

    World.add(world, body);
    ligands.push(body);
    return body;
}

// =============================================================================
// FÍSICA: MOVIMENTO BROWNIANO — O(n)
// =============================================================================

/**
 * Aplica forças brownianas a todos os ligantes.
 *
 * Simula difusão térmica no citosol. A força é aleatória em direção
 * e magnitude, com amortecimento via frictionAir (viscosidade).
 *
 * Complexidade: O(n) — uma operação por ligante por frame.
 */
function applyBrownianMotion() {
    for (const lig of ligands) {
        if (lig.fluidData.phase === 'docked') continue;

        const fx = (Math.random() - 0.5) * CONFIG.BROWNIAN_FORCE;
        const fy = (Math.random() - 0.5) * CONFIG.BROWNIAN_FORCE;
        Body.applyForce(lig, lig.position, { x: fx, y: fy });
    }

    // Proteínas também flutuam suavemente (muito menos)
    for (const prot of proteins) {
        const fx = (Math.random() - 0.5) * CONFIG.BROWNIAN_FORCE * 0.05;
        const fy = (Math.random() - 0.5) * CONFIG.BROWNIAN_FORCE * 0.05;
        Body.applyForce(prot, prot.position, { x: fx, y: fy });
    }
}

// =============================================================================
// TERMODINÂMICA: ATRAÇÃO POR AFINIDADE — O(m)
// =============================================================================

/**
 * Processa resultados de docking e cria molas de atração.
 *
 * Chamado quando o backend retorna resultados do Vina.
 * Se ΔG < AFFINITY_THRESHOLD, cria Constraint invisível (mola)
 * que puxa magneticamente o ligante para a proteína.
 *
 * Complexidade: O(m) onde m = número de matches.
 */
function processDockingResults(results) {
    for (const r of results) {
        const pairKey = `${r.ligand_id}::${r.protein_id}`;
        if (matchedPairs.has(pairKey)) continue;

        if (r.affinity_kcal < CONFIG.AFFINITY_THRESHOLD) {
            const lig = ligands.find(l => l.fluidData.id === r.ligand_id);
            const prot = proteins.find(p => p.fluidData.id === r.protein_id);

            if (!lig || !prot) continue;

            // Transição: floating → attracting
            lig.fluidData.phase = 'attracting';
            lig.fluidData.targetProtein = prot;
            lig.fluidData.affinity = r.affinity_kcal;
            lig.fluidData.color = CONFIG.MATCH_COLOR;

            // Cria mola invisível (Constraint) — atração gradual
            const spring = Constraint.create({
                bodyA: lig,
                bodyB: prot,
                stiffness: CONFIG.SPRING_STIFFNESS,
                damping: 0.02,
                length: 0,  // Quer distância zero (encaixe)
                render: { visible: false },
            });

            World.add(world, spring);
            constraints.push(spring);
            matchedPairs.add(pairKey);

            dockingResults[r.ligand_id] = r;

            console.log(
                `🧲 Match: ${r.ligand_id} → ${r.protein_id} ` +
                `(ΔG = ${r.affinity_kcal.toFixed(1)} kcal/mol)`
            );
        }
    }
}

/**
 * Verifica proximidade e executa "lock-in" (encaixe visual).
 *
 * Complexidade: O(m) — apenas ligantes em fase 'attracting'.
 */
function checkDockingLock() {
    for (const lig of ligands) {
        if (lig.fluidData.phase !== 'attracting') continue;

        const prot = lig.fluidData.targetProtein;
        if (!prot) continue;

        const dx = prot.position.x - lig.position.x;
        const dy = prot.position.y - lig.position.y;
        const dist = Math.sqrt(dx * dx + dy * dy);

        if (dist < CONFIG.LOCK_DISTANCE) {
            // ENCAIXE! Transição: attracting → docked
            lig.fluidData.phase = 'docked';
            Body.setStatic(lig, true);
            Body.setPosition(lig, {
                x: prot.position.x + (Math.random() - 0.5) * 10,
                y: prot.position.y + (Math.random() - 0.5) * 10,
            });

            // Animação de lock (escala pulsa)
            lig.fluidData.scale = 1.5;
            prot.fluidData.glowIntensity = 1.0;

            console.log(`🔒 Docked: ${lig.fluidData.id} → ${prot.fluidData.id}`);
        }
    }
}

// =============================================================================
// RENDERIZAÇÃO p5.js — O(n)
// =============================================================================

/**
 * Setup do p5.js — chamado uma vez.
 */
function setup() {
    const canvas = createCanvas(CONFIG.WIDTH, CONFIG.HEIGHT);
    canvas.parent('fluid-container');
    initPhysics();
}

/**
 * Draw loop do p5.js — chamado a cada frame (~60fps).
 *
 * Complexidade total por frame:
 *   O(n) browniano + O(n) broadphase + O(m) atração + O(n) render = O(n)
 */
function draw() {
    // Atualiza física
    Engine.update(engine, 1000 / 60);
    applyBrownianMotion();
    checkDockingLock();

    // Fundo com trail (efeito fluido)
    background(10, 14, 26, 30);

    // Renderiza proteínas
    for (const prot of proteins) {
        push();
        translate(prot.position.x, prot.position.y);
        rotate(prot.angle);

        const d = prot.fluidData;

        // Glow effect para proteínas com docking ativo
        if (d.glowIntensity > 0) {
            noStroke();
            fill(red(CONFIG.MATCH_COLOR), green(CONFIG.MATCH_COLOR),
                 blue(CONFIG.MATCH_COLOR), d.glowIntensity * 60);
            ellipse(0, 0, d.size * 3, d.size * 3);
            d.glowIntensity *= 0.98;  // Decay
        }

        // Corpo da proteína
        stroke(d.color);
        strokeWeight(2);
        fill(red(d.color), green(d.color), blue(d.color), 40);

        // Desenha vértices
        const verts = prot.vertices;
        beginShape();
        for (const v of verts) {
            vertex(v.x - prot.position.x, v.y - prot.position.y);
        }
        endShape(CLOSE);

        // Label
        noStroke();
        fill(255, 255, 255, 200);
        textAlign(CENTER, CENTER);
        textSize(CONFIG.FONT_SIZE);
        text(d.displayLabel, 0, 0);

        pop();
    }

    // Renderiza ligantes
    for (const lig of ligands) {
        push();
        translate(lig.position.x, lig.position.y);

        const d = lig.fluidData;
        const r = d.radius * d.scale;

        // Animação de escala (decay para normal)
        if (d.scale > 1.0) d.scale *= 0.96;
        if (d.scale < 1.01) d.scale = 1.0;

        // Cor baseada na fase
        const c = color(d.color);
        const alpha = d.phase === 'docked' ? 200 : 255 * d.opacity;

        // Trail de atração
        if (d.phase === 'attracting' && d.targetProtein) {
            stroke(CONFIG.MATCH_COLOR + '40');
            strokeWeight(1);
            line(0, 0,
                d.targetProtein.position.x - lig.position.x,
                d.targetProtein.position.y - lig.position.y);
        }

        // Corpo do ligante
        noStroke();
        fill(red(c), green(c), blue(c), alpha);
        ellipse(0, 0, r * 2, r * 2);

        // Borda luminosa
        noFill();
        stroke(red(c), green(c), blue(c), alpha * 0.5);
        strokeWeight(1);
        ellipse(0, 0, r * 2.3, r * 2.3);

        // Label
        noStroke();
        fill(255, 255, 255, 180);
        textAlign(CENTER, CENTER);
        textSize(CONFIG.FONT_SIZE - 2);
        text(d.displayLabel, 0, r + 12);

        // Badge de afinidade para matches
        if (d.phase === 'docked' || d.phase === 'attracting') {
            const aff = dockingResults[d.id];
            if (aff) {
                fill(CONFIG.MATCH_COLOR);
                textSize(CONFIG.FONT_SIZE - 1);
                text(`${aff.affinity_kcal.toFixed(1)}`, 0, -(r + 10));
            }
        }

        pop();
    }

    // HUD
    drawHUD();
}

/**
 * Desenha o heads-up display com métricas.
 */
function drawHUD() {
    noStroke();
    fill(255, 255, 255, 150);
    textAlign(LEFT, TOP);
    textSize(11);
    text(`Proteínas: ${proteins.length}`, 10, 10);
    text(`Ligantes: ${ligands.length}`, 10, 24);
    text(`Matches: ${matchedPairs.size}`, 10, 38);
    text(`FPS: ${Math.round(frameRate())}`, 10, 52);

    // Legenda
    textAlign(RIGHT, TOP);
    fill(CONFIG.PROTEIN_COLOR);
    text('● Proteína', CONFIG.WIDTH - 10, 10);
    fill(CONFIG.LIGAND_COLORS[0]);
    text('● Ligante', CONFIG.WIDTH - 10, 24);
    fill(CONFIG.MATCH_COLOR);
    text('● Match (ΔG < -7.0)', CONFIG.WIDTH - 10, 38);
}

// =============================================================================
// API: Interface com o Dashboard Dash
// =============================================================================

/**
 * Chamado pelo Dash callback para alimentar dados na simulação.
 * Recebe JSON serializado do backend Python.
 */
window.fluidAPI = {
    addProtein: (id, label, x, y, size) => createProtein(id, label, x, y, size || 50),
    addLigand: (id, label, svg, x, y, r) => createLigand(id, label, svg, x, y, r),
    submitResults: (results) => processDockingResults(results),
    getState: () => ({
        proteins: proteins.length,
        ligands: ligands.length,
        matches: matchedPairs.size,
        fps: Math.round(frameRate ? frameRate() : 0),
    }),
    reset: () => {
        World.clear(world, false);
        proteins = []; ligands = []; constraints = [];
        matchedPairs.clear(); dockingResults = {};
        initPhysics();
    },
};

(() => {
  'use strict';
  // the structure payload is embedded gzip + base64 (protein_viewer.gzip_b64) to keep the single-file report small;
  // it is inflated the first time the explorer starts, so page load and the observers below stay synchronous
  let proteins = null, unpacking = false;
  async function unpack() {
    const packed = Uint8Array.from(atob(document.getElementById('protein-data').textContent.trim()), c => c.charCodeAt(0));
    const inflated = new Blob([packed]).stream().pipeThrough(new DecompressionStream('gzip'));
    proteins = JSON.parse(await new Response(inflated).text()).proteins;
  }
  const el = id => document.getElementById(id);
  const status = text => {el('protein-status').textContent = text;};
  let viewer, protein, pocket, confidence = false, busy = false, surfaces = [], slab;
  const color = a => a.b >= 90 ? '#0053d6' : a.b >= 70 ? '#65cbf3' : a.b >= 50 ? '#ffdb13' : '#ff7d45';
  const controls = ['pv-protein','pv-pocket','pv-whole','pv-focus','pv-confidence'];
  function pocketZoom() {
    const bounds = el('protein-canvas').getBoundingClientRect();
    return 0.35 * Math.min(1, bounds.width / bounds.height);
  }
  function lock(value) {busy = value; controls.forEach(id => {el(id).disabled = value;});}
  function detail() {
    el('pv-title').textContent = protein.gene + (pocket ? ` · Pocket ${pocket.rank}` : ' · No detected pocket');
    const target = el('pv-description'); target.replaceChildren();
    const line = text => {const p = document.createElement('p'); p.textContent = text; target.appendChild(p); return p;};
    if (protein.experimental) {
      el('pv-title').textContent = `${protein.gene} · ${protein.pdb_id}`;
      line(`${protein.gene} + ${protein.ligase} + ${protein.compound}`);
      const p = line('Experimental reference · ');
      const link = document.createElement('a'); link.href=protein.source_url; link.textContent='Primary source'; link.target='_blank'; link.rel='noopener'; p.appendChild(link);
      line(protein.region); line(protein.numbering);
      line(`Annotated residues: ${Object.entries(protein.key_residues).map(([n,a])=>a+n).join(', ') || 'not assigned'}`);
      line('Experimental B-factors are not pLDDT; model-confidence coloring is disabled. Adaptor DDB1 is shown as a faint ribbon.');
      line(protein.scope);
      return;
    }
    const source = line(`${protein.model} · version ${protein.version} · `);
    const link = document.createElement('a'); link.textContent = 'AlphaFold source'; link.href = protein.source_url; link.target = '_blank'; link.rel = 'noopener'; source.appendChild(link);
    line(protein.evidence_type);
    if (pocket) {
      line(`Volume ${pocket.volume_a3} Å³ · Depth ${pocket.max_depth_a} Å · Local mean pLDDT ${pocket.mean_plddt}`);
      line(`Chain ${protein.chain} · residues: ${pocket.residues.join(', ')}`).className = 'pv-residues';
    }
    line(`Detection retained residues with pLDDT ≥ ${protein.detection_min_plddt}. Confidence describes coordinates, not drug binding.`);
  }
  function selectRegion() {
    pocket = protein.experimental
      ? (Object.keys(protein.key_residues).length ? {residues:Object.keys(protein.key_residues).map(Number),kind:'degron'} : null)
      : protein.pockets[Number(el('pv-pocket').value)] || null;
    detail();
  }
  function surfaceStyle() {
    const selected = new Set(pocket ? pocket.residues : []);
    return {opacity:Number(el('pv-opacity').value)/100,
      colorfunc:a => protein.experimental && a.chain === protein.ligase_chain ? '#72a3ce'
        : confidence ? color(a) : a.chain === protein.chain && selected.has(a.resi) ? '#e05454' : '#bcc3cc'};
  }
  async function draw(reset = false) {
    if (busy) return;
    lock(true); status('Building molecular surface…');
    try {
      viewer.removeAllSurfaces(); viewer.setStyle({}, {});
      if (confidence && pocket) viewer.addStyle({chain:protein.chain,resi:pocket.residues},{stick:{color:'#e05454',radius:0.19}});
      surfaces = [];
      if (protein.experimental) {
        const proteinAtoms={chain:[protein.target_chain,protein.ligase_chain],hetflag:false};
        viewer.setStyle({chain:protein.adaptor_chain},{cartoon:{color:'#d0d5db',opacity:.15}});
        viewer.setStyle({resn:protein.ligand_resname},{stick:{color:'#dba737',radius:.3}});
        surfaces.push(await viewer.addSurface($3Dmol.SurfaceType.MS,surfaceStyle(),proteinAtoms,proteinAtoms));
      } else {
        surfaces.push(await viewer.addSurface($3Dmol.SurfaceType.MS, surfaceStyle(), {}, {}));
      }
      if (reset) {
        viewer.zoomTo(pocket ? {chain:protein.chain,resi:pocket.residues} : {});
        viewer.zoom(pocket ? pocketZoom() : 0.85); slab = viewer.getSlab(); el('pv-cut').value = 0;
      }
      viewer.render();
      el('pv-legend').innerHTML = protein.experimental
        ? '<span>Grey: substrate · Blue: CRBN · Gold: compound · Red: annotated degron residues</span>'
        : confidence
        ? '<span>pLDDT: <b style="color:#0053d6">≥90</b> · <b style="color:#269ac6">70–90</b> · <b style="color:#ad9000">50–70</b> · <b style="color:#d95e23">&lt;50</b></span><span>Red sticks: selected pocket residues</span>'
        : '<span><i class="pv-dot" style="background:#bcc3cc"></i>Protein surface</span><span><i class="pv-dot" style="background:#e05454"></i>Selected pocket walls</span>';
      status('Ready · Drag to rotate · Scroll to zoom · Molecular surface calculated from the embedded coordinates.');
    } catch(error) {status('Unable to build 3D surface: ' + error.message);}
    finally {lock(false); el('pv-focus').disabled = !pocket; el('pv-confidence').disabled=!!protein.experimental;}
  }
  async function load() {
    if (busy) return;
    protein = proteins[Number(el('pv-protein').value)];
    if(protein.experimental) {confidence=false; el('pv-confidence').setAttribute('aria-pressed','false');}
    el('pv-focus').textContent=protein.experimental ? 'Focus annotated region' : 'Focus pocket';
    viewer.removeAllModels(); viewer.addModel(protein.pdb,'pdb');
    el('pv-pocket').replaceChildren();
    protein.pockets.forEach((p,i) => el('pv-pocket').add(new Option(`Pocket ${p.rank} · ${p.volume_a3} Å³`,i)));
    if (!protein.pockets.length) el('pv-pocket').add(new Option(protein.experimental ? 'Experimental ternary reference' : 'No pocket detected',''));
    selectRegion(); await draw(true);
  }
  function start() {
    if (viewer || !el('protein-canvas').getBoundingClientRect().width) return;
    if (proteins === null) {
      if (!unpacking) {
        unpacking = true;
        status('Loading the embedded structures...');
        unpack().then(start, error => status('Could not read the embedded structures: ' + error.message));
      }
      return;
    }
    if (!proteins.length) {status('No local protein structures are available.'); return;}
    try {
      // The vendored library creates its surface worker from an in-memory Blob.
      // Coordinates and computation stay local; no CDN or remote worker is used.
      $3Dmol.setSyncSurface(false);
      viewer = $3Dmol.createViewer(el('protein-canvas'),{backgroundColor:'#f3f5f7',antialias:true});
      viewer.setProjection('orthographic');
      proteins.forEach((p,i) => el('pv-protein').add(new Option(p.experimental ? `${p.gene} · ${p.pdb_id} reference` : p.gene,i)));
      el('pv-protein').onchange = load;
      el('pv-pocket').onchange = () => {selectRegion(); draw();};
      el('pv-confidence').onclick = () => {confidence = !confidence; el('pv-confidence').setAttribute('aria-pressed',String(confidence)); draw();};
      el('pv-whole').onclick = () => {viewer.zoomTo(); viewer.zoom(.85); slab=viewer.getSlab(); el('pv-cut').value=0; viewer.render();};
      el('pv-focus').onclick = () => {if(pocket) {viewer.zoomTo({chain:protein.chain,resi:pocket.residues}); viewer.zoom(pocketZoom()); slab=viewer.getSlab(); el('pv-cut').value=0; viewer.render();}};
      el('pv-opacity').oninput = () => {if(!busy) {surfaces.forEach(s=>viewer.setSurfaceMaterialStyle(s,surfaceStyle())); viewer.render();}};
      el('pv-cut').oninput = () => {if(slab && !busy) {const t=Number(el('pv-cut').value)/100; viewer.setSlab(slab.near + (slab.far-slab.near)*t*.9,slab.far); viewer.render();}};
      new ResizeObserver(() => {if(el('protein-canvas').getBoundingClientRect().width) viewer.resize();}).observe(el('protein-canvas'));
      load();
    } catch(error) {status('3D needs WebGL in your browser. The tables remain available. ' + error.message);}
  }
  new IntersectionObserver(entries => {if(entries.some(e => e.isIntersecting)) start();}).observe(el('protein-canvas'));
})();

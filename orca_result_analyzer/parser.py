class OrcaParser:
    def __init__(self):
        self.filename = ""
        self.raw_content = ""
        self.lines = []
        self.data = {
            "converged": False,
            "scf_energy": None,
            "atoms": [],
            "coords": [],
            "charge": 0,
            "mult": 1,
            "frequencies": [], # List of dicts: {freq, ir, raman, vector}
            "excitation_energies": [], # TDDFT
            "dipoles": [],
            "nmr_shielding": [],
            "charges": {} # Type -> List
        }

    def load_from_memory(self, content, filename=""):
        self.filename = filename
        self.raw_content = content
        self.lines = content.splitlines()
        self.parse_all()

    def parse_all(self):
        self.parse_basic()
        self.parse_trajectory()
        self.parse_frequencies()
        self.parse_thermal()
        self.parse_mo()
        self.parse_mo_coeffs()
        self.parse_charges()
        self.parse_dipole()
        self.parse_tddft()
        self.parse_gradients()
        self.parse_nmr()
        self.parse_basis_set()
        
    def parse_basic(self):
        """Parse basic info: SCF Energy, Convergence, Geometry."""
        for i, line in enumerate(self.lines):
            line = line.strip()
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    self.data["scf_energy"] = float(line.split()[-1])
                except: pass
            if "Total Charge" in line and "Multiplicity" in line:
                try:
                     parts = line.split()
                     self.data["charge"] = int(parts[-4])
                     self.data["mult"] = int(parts[-1])
                except: pass
            if "SCF CONVERGED AFTER" in line or "Optimization converged" in line or "HURRAY" in line:
                self.data["converged"] = True
                
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                # Read geometry
                self.data["atoms"] = []
                self.data["coords"] = []
                curr = i + 2
                while curr < len(self.lines):
                    l_geo = self.lines[curr].strip()
                    if not l_geo or "-------" in l_geo: break
                    parts = l_geo.split()
                    if len(parts) >= 4:
                        self.data["atoms"].append(parts[0])
                        self.data["coords"].append([float(parts[1]), float(parts[2]), float(parts[3])])
                    curr += 1

    def parse_trajectory(self):
        """Parse Optimization and Scan Trajectories."""
        self.data["scan_steps"] = []
        
        # Look for "RELAXED SURFACE SCAN STEP" or "GEOMETRY OPTIMIZATION CYCLE"
        # And capture Energy + Geometry
        # Actually usually Step header -> Energy -> ... -> Coordinates
        
        current_step = {}
        in_step = False
        
        # Helper to find coords after a header
        def read_coords_from(idx):
            atoms = []
            coords = []
            curr = idx + 2 # Skip header and rule
            # ORCA output for coords in opt steps usually:
            # "CARTESIAN COORDINATES (ANGSTROEM)"
            # But during Opt it prints "CARTESIAN COORDINATES (A.U.)" or similar?
            # Actually standard output prints "CARTESIAN COORDINATES (ANGSTROEM)" often.
            
            # Let's search forward for coordinates
            limit = 500 # Search limit
            found_coords = False
            for k in range(limit):
                if idx + k >= len(self.lines): break
                line = self.lines[idx+k].strip()
                if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                    c_idx = idx + k + 2
                    found_coords = True
                    while c_idx < len(self.lines):
                         cl = self.lines[c_idx].strip()
                         if not cl or "-------" in cl: break
                         parts = cl.split()
                         if len(parts) >= 4:
                             atoms.append(parts[0])
                             coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                         c_idx += 1
                    break
            return atoms, coords, found_coords

        for i, line in enumerate(self.lines):
            # Scan Step Header
            if "RELAXED SURFACE SCAN STEP" in line:
                step_idx = 0
                try: step_idx = int(line.split()[-1])
                except: pass
                
                # Now find Energy and Coords associated with this step
                # Usually prints "FINAL SINGLE POINT ENERGY" shortly after
                # But that might be reused by parse_basic (last one).
                # We need the one IN this block.
                
                # Let's just grab the NEXT Energy and Coords
                # Find energy
                en = 0.0
                for k in range(i, len(self.lines)):
                    if "RELAXED SURFACE SCAN STEP" in self.lines[k] and k > i: break # Next step
                    if "FINAL SINGLE POINT ENERGY" in self.lines[k]:
                         try: en = float(self.lines[k].split()[-1])
                         except: pass
                         break
                
                atoms, coords, found = read_coords_from(i)
                if found:
                    self.data["scan_steps"].append({
                        'step': step_idx,
                        'energy': en,
                        'atoms': atoms,
                        'coords': coords
                    })
                    
            elif "GEOMETRY OPTIMIZATION CYCLE" in line:
                cycle_idx = 0
                try: cycle_idx = int(line.split()[-2]) # Cycle 1
                except: pass
                
                # Find Energy
                en = 0.0
                for k in range(i, len(self.lines)):
                    if "GEOMETRY OPTIMIZATION CYCLE" in self.lines[k] and k > i: break
                    if "FINAL SINGLE POINT ENERGY" in self.lines[k]:
                         try: en = float(self.lines[k].split()[-1])
                         except: pass
                         break
                         
                atoms, coords, found = read_coords_from(i)
                if found:
                     self.data["scan_steps"].append({
                        'step': cycle_idx,
                        'energy': en,
                        'atoms': atoms,
                        'coords': coords
                    })
        
    def parse_mo_coeffs(self):
        self.data["mo_coeffs"] = {} # mo_idx -> { 'coeffs': list, 'energy': float, 'occ': float, 'spin': 'alpha'/'beta' }
        
        # Blocks to look for:
        # "MOLECULAR ORBITALS" (RHF)
        # "SPIN UP ORBITALS" (UHF Alpha)
        # "SPIN DOWN ORBITALS" (UHF Beta)
        
        current_spin = "alpha" # Default
        
        start_indices = []
        for i, line in enumerate(self.lines):
            if "MOLECULAR ORBITALS" in line and "------------------" in self.lines[i+1]:
                start_indices.append((i, "alpha"))
            elif "SPIN UP ORBITALS" in line and "----------------" in self.lines[i+1]:
                start_indices.append((i, "alpha"))
            elif "SPIN DOWN ORBITALS" in line and "------------------" in self.lines[i+1]:
                start_indices.append((i, "beta"))
                
        if not start_indices: return
        
        for start_idx, spin in start_indices:
            curr = start_idx + 2
            
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line: 
                    curr += 1
                    continue
                if "TIMINGS" in line: break 
                if "--------" in line: 
                    curr += 1
                    continue
                if "ORBITALS" in line and "--------" in self.lines[curr+1]: break # Next block
                
                parts = line.split()
                if not parts: 
                    curr+=1
                    continue
                    
                # Header check: Integers? "0   1   2..."
                is_header = False
                try:
                    # Check first few items
                    if all(p.isdigit() for p in parts):
                         indices = [int(p) for p in parts]
                         is_header = True
                except: pass
                
                if is_header:
                    current_mos = [int(p) for p in parts]
                    
                    # Init storage dicts
                    for idx in current_mos:
                        if idx not in self.data["mo_coeffs"]:
                             self.data["mo_coeffs"][idx] = {'coeffs': [], 'spin': spin, 'energy': 0.0, 'occ': 0.0}
                    
                    # Try to parse Energy / Occ lines immediately following
                    # Usually:
                    #                  -10.000    -5.000
                    #                   2.0000     1.000
                    try:
                        next_lines = [self.lines[curr+1].strip(), self.lines[curr+2].strip()]
                        # Check if numbers
                        vals1 = next_lines[0].split()
                        vals2 = next_lines[1].split()
                        
                        have_energy_occ = False
                        if len(vals1) == len(current_mos) and len(vals2) == len(current_mos):
                             # Assume line 1 is Energy (Eh), Line 2 is Occ
                             # Verify they look like floats
                             try:
                                 energies = [float(v) for v in vals1]
                                 occs = [float(v) for v in vals2]
                                 for k, idx in enumerate(current_mos):
                                     self.data["mo_coeffs"][idx]['energy'] = energies[k]
                                     self.data["mo_coeffs"][idx]['occ'] = occs[k]
                                 have_energy_occ = True
                                 curr += 2 # Skip these 2 lines
                             except: pass
                    except: pass
                             
                    curr += 1
                    continue
                
                # Coefficient line: "0   C  1s   0.000 ..."
                if len(parts) >= 3:
                     # Check if start is int index (0) and sym (C)
                     try:
                         atom_idx = int(parts[0])
                         sym = parts[1]
                         orb = parts[2]
                         val_strs = parts[3:]
                         
                         # Sometimes sym is attached? 0C?
                         # Usually ORCA puts spaces: "   0   C   1s"
                         
                         # Check if remaining parts match number of current MOs
                         if len(val_strs) == len(current_mos):
                              for k, v_str in enumerate(val_strs):
                                  mo_idx = current_mos[k]
                                  val = float(v_str)
                                  if abs(val) > 0.01:
                                      self.data["mo_coeffs"][mo_idx]['coeffs'].append({
                                          "atom_idx": atom_idx,
                                          "sym": sym,
                                          "orb": orb,
                                          "coeff": val
                                      })
                     except: pass
                curr += 1
        
    def parse_charges(self):
        self.data["charges"] = {} 
        
        # Look for Mulliken, Loewdin, Hirshfeld
        # MULLIKEN ATOMIC CHARGES
        # LOEWDIN ATOMIC CHARGES
        # HIRSHFELD ANALYSIS
        
        blocks = {
            "Mulliken": "MULLIKEN ATOMIC CHARGES",
            "Loewdin": "LOEWDIN ATOMIC CHARGES",
            "Hirshfeld": "HIRSHFELD ANALYSIS"
        }
        
        for name, header in blocks.items():
            for i, line in enumerate(self.lines):
                if header in line:
                    # Found block
                    # Skip until "-----------------"
                    # Then read
                    curr = i + 1
                    while curr < len(self.lines):
                         if "----------------" in self.lines[curr]:
                             curr += 1
                             break
                         curr += 1
                    
                    found_charges = []
                    while curr < len(self.lines):
                        l = self.lines[curr].strip()
                        if not l: 
                            if found_charges: break # Done reading
                            curr += 1
                            continue
                        if "Sum of atomic charges" in l: break
                        if "----------------" in l: break
                        
                        parts = l.split()
                        # Format: idx Atom Charge ...
                        # 0 C : -0.123
                        # OR 0 C -0.123
                        if len(parts) >= 3:
                             # Check if first is int
                             try:
                                 idx = int(parts[0].strip(':'))
                                 sym = parts[1]
                                 val = float(parts[-1]) # Charge is usually last or 3rd
                                 # In Hirshfeld: "  0   C     -0.033606     0.054394" (Charge, Spin)
                                 # In Mulliken: "  0 C :   -0.123 "
                                 
                                 # Refine Parsing based on name matches
                                 if name == "Hirshfeld":
                                     # 0  C   Charge   Spin
                                     val = float(parts[2])
                                 else:
                                     # Mulliken/Loewdin
                                     # 0 C : -0.123
                                     # Check if parts[2] is ':'
                                     if parts[2] == ':':
                                         val = float(parts[3])
                                     else:
                                         val = float(parts[2])
                                         
                                 found_charges.append({'atom_idx': idx, 'atom_sym': sym, 'charge': val})
                             except: pass
                        curr += 1
                        
                    if found_charges:
                        self.data["charges"][name] = found_charges

    def parse_scan(self):
        """Alias for parse_trajectory."""
        self.parse_trajectory()

    def parse_gradient(self):
        """Alias for parse_gradients."""
        self.parse_gradients()

    def parse_gradients(self):
        """Parse Cartesian Gradients"""
        self.data["gradients"] = []
        found = False
        start_idx = -1
        for i, line in enumerate(self.lines):
            if "CARTESIAN GRADIENTS" in line and "-------" not in line:
                start_idx = i
                found = True
        
        if not found: return
        
        curr = start_idx + 2
        while curr < len(self.lines):
            line = self.lines[curr].strip()
            if "-------" in line: break
            if not line: break
            
            parts = line.split()
            # 1 C : X Y Z
            if len(parts) >= 6:
                try:
                    idx = int(parts[0])
                    sym = parts[1]
                    vx = float(parts[3])
                    vy = float(parts[4])
                    vz = float(parts[5])
                    self.data["gradients"].append({
                        'atom_idx': idx, 
                        'atom_sym': sym, 
                        'vector': [vx, vy, vz]
                    })
                except: pass
            
            curr += 1

    def parse_dipole(self):
        # Look for "TOTAL DIPOLE MOMENT"
        # -------------
        # DIPOLE MOMENT
        # -------------
        # ...
        # Total Dipole Moment    :   1.234   2.345   3.456
        # Magnitude (Debye)      :   4.567
        
        self.data["dipole"] = None
        
        start_idx = -1
        # Reverse search might be safer for final dipole
        # But usually explicitly labeled "DIPOLE MOMENT" block.
        # Let's search from end or find all and take last.
        
        candidates = []
        for i, line in enumerate(self.lines):
            if "Total Dipole Moment" in line and ":" in line:
                 candidates.append(i)
                 
        if not candidates: return
        
        # Take the last occurrence
        idx = candidates[-1]
        line = self.lines[idx]
        # Format: Total Dipole Moment    :    -0.00000     0.00000     0.00000
        parts = line.split(":")
        if len(parts) > 1:
            try:
                vec_str = parts[1].strip().split()
                if len(vec_str) >= 3:
                     x, y, z = float(vec_str[0]), float(vec_str[1]), float(vec_str[2])
                     
                     # Look for magnitude
                     # Next line maybe?
                     # Magnitude (Debye)      :      0.00000
                     mag = 0.0
                     if idx + 1 < len(self.lines):
                         line2 = self.lines[idx+1]
                         if "Magnitude" in line2 and ":" in line2:
                             mag = float(line2.split(":")[1].strip())
                         else:
                             # Calculate manually
                             import math
                             mag = math.sqrt(x*x + y*y + z*z)
                             
                     self.data["dipole"] = {
                         "vector": (x, y, z),
                         "magnitude": mag
                     }
            except: pass

        self.parse_nmr()
        self.parse_charges()
        self.parse_charges()
        self.parse_scan()
        self.parse_gradient()
        self.parse_basis_set()
        
    def parse_basis_set(self):
        self.data["basis_set_shells"] = [] # List of {type, center, exps, coeffs}
        
        # 1. Parse GTO Definitions
        # Look for "BASIS SET IN INPUT FORMAT"
        start_idx = -1
        for i, line in enumerate(self.lines):
            if "BASIS SET IN INPUT FORMAT" in line:
                start_idx = i
                break
        
        if start_idx == -1: return
        
        gto_defs = {} # Sym -> List of shelldefs {type: 'S', exps: [], coeffs: []}
        
        curr = start_idx + 2
        current_sym = None
        current_shells = []
        
        # Helper to map type string to int
        type_map = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4} 
        
        while curr < len(self.lines):
            line = self.lines[curr].strip()
            if "Basis set for element" in line:
                 # "Basis set for element : C"
                 pass # Just info, usually followed by NewGTO
            
            if line.startswith("NewGTO"):
                parts = line.split()
                if len(parts) >= 2:
                    current_sym = parts[1]
                    current_shells = []
            
            elif line.startswith("EndGTO"):
                if current_sym:
                    gto_defs[current_sym] = current_shells
                    current_sym = None
                    
            elif line.startswith("-----"): # End of block
                 break
            
            elif current_sym is not None:
                # Inside a GTO block
                # Format:
                # S   3
                #  1.23   0.45
                #  ...
                parts = line.split()
                if len(parts) >= 2 and parts[0] in type_map:
                    # New shell header: Type  NPrims
                    stype_str = parts[0]
                    n_prim = int(parts[1])
                    
                    stype = type_map[stype_str]
                    exps = []
                    coeffs = []
                    
                    # Read primitives
                    for _ in range(n_prim):
                        curr += 1
                        if curr >= len(self.lines): break
                        pline = self.lines[curr].strip()
                        pparts = pline.split()
                        # Often: Index Exp Coeff (ORCA sometimes adds index 1, 2...)
                        # Or just: Exp Coeff
                        # Check format.
                        # ORCA input format: "  1   123.45   0.123"
                        try:
                            if len(pparts) == 3:
                                exps.append(float(pparts[1]))
                                coeffs.append(float(pparts[2]))
                            elif len(pparts) == 2:
                                exps.append(float(pparts[0]))
                                coeffs.append(float(pparts[1]))
                        except: pass
                    
                    current_shells.append({
                        "type": stype,
                        "exps": exps,
                        "coeffs": coeffs
                    })
            
            curr += 1

        # 2. Expand to full molecule (Atom Centers)
        # We need self.data["atoms"] and self.data["coords"]
        if not self.data["atoms"] or not self.data["coords"]:
             # Maybe geometry parsing failed or hasn't run? 
             # It should have run in parse_basic -> parse_all sequence.
             return

        # Unit conversion: Parser stores coords in Angstrom.
        # MO Engine expects Bohr (usually).
        # Let's verify what `mo_engine` does. 
        # `mo_engine` says "atoms_coords: list/array of (x,y,z) in Angstrom" for CalcWorker init,
        # BUT `BasisSetEngine` expects "center: [x, y, z] (BOHR)".
        # CalcWorker does conversion: coords_bohr = self.atoms_coords * ANG_TO_BOHR.
        # But CalcWorker does NOT pass `shells` to Engine. It passes `mo_coeffs` and generic atoms.
        # WAIT. `mo_engine` needs `shells` initialized in `BasisSetEngine`.
        # `CalcWorker` takes `engine` as argument.
        # So WE need to initialize `BasisSetEngine` with shelsl in BOHR.
        
        BOHR_TO_ANG = 0.529177249
        ANG_TO_BOHR = 1.0 / BOHR_TO_ANG
        
        full_shells = []
        
        for idx, sym in enumerate(self.data["atoms"]):
            # Look up def for this symbol
            # ORCA might have specific basis for specific atoms "H:1"?
            # For now assume mostly element based.
            
            defs = gto_defs.get(sym)
            if not defs:
                # Try checking if symbol has numbers? "0C"? 
                # Parser usually cleans sym.
                continue
                
            coord_ang = self.data["coords"][idx]
            coord_bohr = [c * ANG_TO_BOHR for c in coord_ang]
            
            import numpy as np
            
            for sh_def in defs:
                full_shells.append({
                    "type": sh_def["type"],
                    "center": np.array(coord_bohr),
                    "exps": np.array(sh_def["exps"]),
                    "coeffs": np.array(sh_def["coeffs"])
                })
        
        self.data["basis_set_shells"] = full_shells

        self.data["gradients"] = [] # List of {atom_idx, atom_sym, vector: (x,y,z)}
        
        # Look for CARTESIAN GRADIENT
        # ------------------
        # CARTESIAN GRADIENT
        # ------------------
        #
        #   1   C   :   -0.000003588    0.000001235    0.000000985
        
        start_idx = -1
        for i, line in enumerate(self.lines):
            if "CARTESIAN GRADIENT" in line:
                start_idx = i
                # Could be multiple (optimization steps). Use last one?
                # Usually scan implies last one.
        
        if start_idx != -1:
             curr = start_idx + 3 # Skip header lines (adjusted during run if needed)
             
             while curr < len(self.lines):
                 line = self.lines[curr].strip()
                 if "Difference to the" in line or "Total Gradient" in line: break # End of block
                 if not line: 
                     curr += 1
                     continue
                 
                 parts = line.split()
                 # Format: Index Sym : X Y Z
                 # 0 C : ...
                 # ORCA sometimes:   0   C   : ...
                 if len(parts) >= 6 and parts[2] == ":":
                     try:
                         idx = int(parts[0])
                         sym = parts[1]
                         gx = float(parts[3])
                         gy = float(parts[4])
                         gz = float(parts[5])
                         self.data["gradients"].append({
                             "atom_idx": idx,
                             "atom_sym": sym,
                             "vector": (gx, gy, gz)
                         })
                     except: pass
                 elif len(parts) >= 5: # Maybe different format
                      # Try parsing last 3 as floats
                      try:
                          gx, gy, gz = float(parts[-3]), float(parts[-2]), float(parts[-1])
                          sym = parts[1]
                          idx = int(parts[0])
                          self.data["gradients"].append({
                             "atom_idx": idx,
                             "atom_sym": sym,
                             "vector": (gx, gy, gz)
                         })
                      except: pass
                 
                 curr += 1

        self.data["scan_steps"] = [] # List of {step: int, energy: float, atoms: [], coords: []}
        
        # Two cases: 
        # 1. "Relaxed Surface Scan" where ORCA prints "The current geometry is:" after each step.
        # 2. Sequential optimization steps (e.g. trajectory) where getting every "CARTESIAN COORDINATES" block helps.
        
        # Strategy:
        # Find "Relaxed Surface Scan" header to confirm it IS a scan?
        # Or just collect all Geometries + Energies if multiple exist?
        
        # Let's try to collect all "CARTESIAN COORDINATES (ANGSTROEM)" blocks.
        # But we need corresponding Energy.
        
        # ORCA "Relaxed Surface Scan" often prints:
        # "The current geometry is:" 
        # followed by XYZ format.
        
        is_scan = False
        for line in self.lines:
            if "Relaxed Surface Scan" in line:
                is_scan = True
                break
                
        # Even if not explicit scan, maybe trajectory?
        # Let's look for "CARTESIAN COORDINATES (ANGSTROEM)" and "FINAL SINGLE POINT ENERGY"
        # The order matters. Usually Geometry then Energy or vice versa in output stream.
        
        # New approach:
        # Iterate and collect blocks.
        
        steps = []
        
        current_coords = []
        current_atoms = []
        
        # Regex or simple search?
        # ORCA 5:
        # ---------------------------------
        # CARTESIAN COORDINATES (ANGSTROEM)
        # ---------------------------------
        # ... atoms ...
        #
        # ...
        # FINAL SINGLE POINT ENERGY ...
        
        # But in a scan, there are multiple.
        
        coord_indices = []
        for i, line in enumerate(self.lines):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                coord_indices.append(i)
                
        # Parse each coordinate block
        candidates = []
        for start_idx in coord_indices:
            atoms = []
            coords = []
            curr = start_idx + 2
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line: 
                    curr += 1
                    continue
                parts = line.split()
                if len(parts) >= 4: # Sym X Y Z
                    try:
                        sym = parts[0]
                        # Check if first part is valid atom or just text
                        if not sym[0].isalpha(): break # Header or separator
                        
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        atoms.append(sym)
                        coords.append((x, y, z))
                    except: break
                else:
                    break
                curr += 1
            
            if atoms:
                 candidates.append({'atoms': atoms, 'coords': coords, 'line': start_idx, 'energy': None})
                 
        # Now try to attach energies.
        # Usually Energy is printed lines AFTER coords or BEFORE?
        # In optimization: Coords printed, then "FINAL SINGLE POINT ENERGY" follows.
        
        # For each candidate, search forward for Energy until next candidate or end.
        
        for k in range(len(candidates)):
            cand = candidates[k]
            start_search = cand['line']
            end_search = candidates[k+1]['line'] if k + 1 < len(candidates) else len(self.lines)
            
            # Search for energy in [start_search, end_search]
            best_energy = None
            
            # ORCA Scan often prints "FINAL SINGLE POINT ENERGY" at end of step.
            # But there might be intermediate SCF energies.
            # We want the "Converged" one usually? 
            # "Optimization Converged" usually precedes the FINAL energy printout in that cycle.
            
            # Search logic: find last "FINAL SINGLE POINT ENERGY" in this block?
            # Or "Total Energy       :"
            
            for r in range(start_search, end_search):
                line = self.lines[r]
                if "FINAL SINGLE POINT ENERGY" in line:
                    try:
                        val = float(line.split()[-1])
                        best_energy = val
                    except: pass
            
            cand['energy'] = best_energy
            
        # Filter valid steps (must have energy? or at least coords)
        # If it's a scan, we expect changing energies.
        
        valid_steps = [c for c in candidates if c['energy'] is not None]
        
        # If too few (just 1), maybe just a single point calc?
        # But user wants Scan Graph.
        if len(valid_steps) > 1:
            for i, s in enumerate(valid_steps):
                self.data["scan_steps"].append({
                    "step": i+1,
                    "energy": s['energy'],
                    "atoms": s['atoms'],
                    "coords": s['coords']
                })

        self.data["charges"] = {} # type -> list of {atom_idx, atom_sym, charge}
        
        # MULLIKEN ATOMIC CHARGES

        mulliken_start = -1
        # LOEWDIN ATOMIC CHARGES
        loewdin_start = -1
        # HIRSHFELD ANALYSIS
        hirshfeld_start = -1
        
        for i, line in enumerate(self.lines):
            if "MULLIKEN ATOMIC CHARGES" in line: mulliken_start = i
            if "LOEWDIN ATOMIC CHARGES" in line: loewdin_start = i
            if "HIRSHFELD ANALYSIS" in line: hirshfeld_start = i
            
        def parse_block(start_idx, header_lines=2):
            res = []
            if start_idx == -1: return res
            curr = start_idx + header_lines
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line: break
                if "Sum of atomic charges" in line: break
                if "----------------" in line: break
                
                parts = line.split()
                # Format: 0 C : -0.123
                # Or: 0 C -0.123
                if len(parts) >= 3:
                     # Usually: Index Symbol : Charge
                     # 0 O : -0.567
                     try:
                         idx_str = parts[0].strip(":")
                         sym = parts[1].strip(":")
                         val_str = parts[-1] 
                         # Careful with ":"
                         if parts[2] == ":": val_str = parts[3]
                         
                         res.append({
                             "atom_idx": int(idx_str),
                             "atom_sym": sym,
                             "charge": float(val_str)
                         })
                     except: pass
                curr += 1
            return res

        self.data["charges"]["Mulliken"] = parse_block(mulliken_start)
        self.data["charges"]["Loewdin"] = parse_block(loewdin_start)
        # Hirshfeld usually has a different format, verify if needed. 
        # Usually: "  Atom     Charge  Spin" loop
        # For now, let's stick to Mulliken/Loewdin as they are most common.

        
    def parse_nmr(self):
        self.data["nmr_shielding"] = [] # List of {atom_idx, atom_sym, shielding, shift=None}
        
        # Look for "CHEMICAL SHIELDING SUMMARY (ppm)" or individual nucleus blocks
        # ORCA output typically has a summary at the end of properties
        
        summary_start = -1
        # Try to find the summary table first
        for i, line in enumerate(self.lines):
            if "CHEMICAL SHIELDING SUMMARY (ppm)" in line:
                summary_start = i
                
        if summary_start != -1:
            curr = summary_start + 5 # Skip headers
            #  N Nucleus  Shielding
            #  0  H       31.456
            #  ...
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line: break
                if "---------------" in line: break
                
                parts = line.split()
                # Format: Index Nucleus Shielding
                if len(parts) >= 3:
                     try:
                         idx = int(parts[0])
                         sym = parts[1]
                         val = float(parts[2])
                         self.data["nmr_shielding"].append({
                             "atom_idx": idx,
                             "atom_sym": sym,
                             "shielding": val
                         })
                     except: pass
                curr += 1
                
        else:
            # Fallback: Parse individual blocks "CHEMICAL SHIELDING TENSOR (ppm)"
            # Nucleus   0C :
            # ...
            # Isotropic   =    123.456
            
            # This is more complex as it's scattered.
            # Let's search for "Nucleus" and "Isotropic" lines if summary not found.
            pass # Summary is almost always there in recent ORCA versions for NMR jobs

        
    def parse_tddft(self):
        self.data["tddft"] = [] # List of {state, energy_ev, energy_nm, osc}
        
        # Look for "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
        # This is the most reliable summary table in ORCA 5.x
        
        # Sample Format:
        #      Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ   
        #                       (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)  
        # ----------------------------------------------------------------------------------------------------
        #   0-1A  ->  2-1A    ... (values) ...
        
        found_block = False
        start_idx = -1
        
        for i, line in enumerate(self.lines):
            if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
                start_idx = i
                found_block = True
                break
                
        if found_block:
            curr = start_idx + 1
            # Skip until we hit data
            # Usually there is a header, then unit line, then separator line
            # Let's just look for the separator line "---------"
            
            while curr < len(self.lines):
                if "----------------" in self.lines[curr]:
                    curr += 1
                    break
                curr += 1
                
            # Now we are at unit line? Or Header?
            # actually usually:
            # Header
            # Units
            # Separator
            # Data
            
            # So if we found one separator, might need to find the NEXT one depending on where we started
            # Let's strictly skip header lines until we find a line starting with "0-" or similar, or just parse based on column count and type
            
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line: 
                    curr += 1
                    continue
                if "----------------" in line: 
                    curr += 1 
                    continue # Skip separators
                    
                if "SOC CORRECTED" in line: break
                if "CD SPECTRUM" in line: break
                
                parts = line.split()
                # Expected: 0-1A -> 2-1A Ev Cm-1 Nm Osc ...
                # 0: 0-1A
                # 1: ->
                # 2: 2-1A (Target State)
                # 3: eV
                # 4: cm-1
                # 5: nm
                # 6: fosc
                
                if len(parts) >= 7 and "->" in parts:
                    try:
                        # Parse State ID from "2-1A" -> 2
                        target_state = parts[2]
                        # remove -1A etc
                        # just take leading digits?
                        state_idx_str = ""
                        for char in target_state:
                            if char.isdigit(): state_idx_str += char
                            else: break
                        
                        state = int(state_idx_str) if state_idx_str else 0
                        
                        ev = float(parts[3])
                        nm = float(parts[5])
                        osc = float(parts[6])
                        
                        self.data["tddft"].append({
                            "state": state,
                            "energy_ev": ev,
                            "energy_nm": nm,
                            "osc": osc
                        })
                    except: pass
                # Fallback for older formats or different tables if "->" not present?
                # Older ORCA might be: 1   26755.5    373.7    0.0000000
                elif len(parts) >= 4 and parts[0].isdigit():
                     try:
                        state = int(parts[0])
                        # Check units. Usually cm-1, nm, osc
                        # But verifying via header is hard.
                        # Heuristic: 2nd col usually large (cm-1), 3rd col ~hundreds (nm), 4th small (osc)
                        v2 = float(parts[2])
                        v3 = float(parts[3])
                        
                        # Just assume standard order: State, Energy(cm-1), Wavelength(nm), Osc
                        nm = v2
                        osc = v3
                        ev = 1239.84193 / nm if nm else 0
                        
                        self.data["tddft"].append({
                             "state": state,
                             "energy_ev": ev,
                             "energy_nm": nm,
                             "osc": osc
                        })
                     except: pass
                     
                curr += 1
                
        # Parse CD Spectrum if available
        # Header: CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
        # Cols: Transition, Energy(eV), Energy(cm-1), Wavelength(nm), R (1e40*cgs), MX, MY, MZ
        
        cd_start = -1
        for i, line in enumerate(self.lines):
            if "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
                cd_start = i
                break
                
        if cd_start != -1:
            curr = cd_start + 1
            while curr < len(self.lines):
                if "----------------" in self.lines[curr]:
                    curr += 1
                    break
                curr += 1
                
            # Skip potential second separator or units
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue
                if "----------------" in line:
                    curr += 1
                    continue
                
                parts = line.split()
                # 0-1A -> 2-1A  eV  cm-1  nm  R ...
                
                if len(parts) >= 6 and "->" in parts:
                    try:
                        # Extract state from "2-1A"
                        target_state = parts[2]
                        state_idx_str = ""
                        for char in target_state:
                            if char.isdigit(): state_idx_str += char
                            else: break
                        state = int(state_idx_str) if state_idx_str else 0
                        
                        # R is usually at index 6 if standard format
                        # Transition(0,1,2), eV(3), cm(4), nm(5), R(6)
                        rot_str = float(parts[6])
                        
                        # Find matching state in self.data["tddft"]
                        for item in self.data["tddft"]:
                            if item["state"] == state:
                                item["rotatory_strength"] = rot_str
                                break
                    except: pass
                curr += 1
        
    def parse_mo(self):
        self.data["mos"] = [] # List of dict: {id, occ, energy, spin}
        
        # Look for "ORBITAL ENERGIES"
        start_line = -1
        for i, line in enumerate(self.lines):
            if "ORBITAL ENERGIES" in line:
                start_line = i
        
        if start_line != -1:
            curr = start_line + 4 # Skip headers
            # Format:
            #   NO   OCC          E(Eh)            E(eV)
            #    0   2.0000      -19.125026      -520.4184
            #   ...
            
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue
                if "MULLIKEN POPULATION ANALYSIS" in line: break
                
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        idx = int(parts[0])
                        occ = float(parts[1])
                        eh = float(parts[2])
                        ev = float(parts[3])
                        
                        # Check spin? ORCA lists alpha/beta separately if UKS.
                        # Usually "SPIN UP ORBITALS" and "SPIN DOWN" blocks if unrestricted.
                        # For now assume restricted or just list all found.
                        
                        self.data["mos"].append({
                            "id": idx,
                            "occ": occ,
                            "energy_eh": eh,
                            "energy_ev": ev
                        })
                    except: pass
                curr += 1

    def parse_thermal(self):
        self.data["thermal"] = {}
        # ORCA Thermochem block
        # Look for "THERMOCHEMISTRY AT 298.15 K"
        
        start_line = -1
        for i, line in enumerate(self.lines):
            if "THERMOCHEMISTRY AT" in line:
                start_line = i
                
        if start_line != -1:
            curr = start_line
            # Scan down for specific keywords
            # "Total thermal energy                  ..."
            # "Total Enthalpy                        ..."
            # "Final Entropy                         ..."
            # "Final Gibbs free energy               ..."
            
            thermal_keys = {
                "Zero point energy": "zpe",
                "Total thermal energy": "thermal_energy",
                "Total Enthalpy": "enthalpy",
                "Final Entropy": "entropy",
                "Final Gibbs free energy": "gibbs"
            }
            
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue
                if "Timings for individual modules" in line: break # Safe stop
                if curr > start_line + 100: break # Safety
                
                for key, val_key in thermal_keys.items():
                    if key in line:
                         parts = line.split()
                         try:
                             # usually last value is Hartrees, maybe "..." before it
                             val = float(parts[-1])
                             if "Entropy" in key:
                                 # Entropy line might be different unit or format?
                                 # ORCA: "Final Entropy                      ...    0.000000 Eh" ?
                                 # Usually printed as:
                                 # "Total Entropy correction           ...   -0.123456 Eh"
                                 pass
                             self.data["thermal"][val_key] = val
                         except: pass
                curr += 1
            
        # 2. Final Energy
        for line in reversed(self.lines):
            target = "FINAL SINGLE POINT ENERGY"
            if target in line:
                try:
                    parts = line.split()
                    val = float(parts[-1])
                    self.data["scf_energy"] = val
                    break
                except: pass

        # 3. Geometry (Last occurrence)
        coord_start = -1
        for i, line in enumerate(self.lines):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                coord_start = i
        
        if coord_start != -1:
            self.data["atoms"] = []
            self.data["coords"] = []
            curr = coord_start + 2
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue
                    
                parts = line.split()
                if len(parts) >= 4:
                    sym = parts[0]
                    try:
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        self.data["atoms"].append(sym)
                        self.data["coords"].append((x, y, z))
                    except ValueError:
                        break
                else:
                    break
                curr += 1

    def parse_frequencies(self):
        self.data["frequencies"] = []
        
        # 1. Frequencies
        freq_start = -1
        for i, line in enumerate(self.lines):
            if "VIBRATIONAL FREQUENCIES" in line:
                freq_start = i
                # Don't break immediately, could be multiple? usually last one matters or first?
                # In optimization + freq, it's at end.
        
        if freq_start != -1:
            curr = freq_start + 4
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "NORMAL MODES" in line: break
                if not line:
                    curr += 1
                    continue
                
                if ":" in line and "cm**-1" in line:
                    parts = line.split()
                    try:
                        val = float(parts[1])
                        self.data["frequencies"].append({"freq": val, "ir": 0.0, "raman": 0.0, "vector": []})
                    except: pass
                curr += 1

        # 2. IR Spectrum
        ir_start = -1
        for i, line in enumerate(self.lines):
            if "IR SPECTRUM" in line:
                ir_start = i
        
        if ir_start != -1:
            curr = ir_start + 5
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "The first frequency" in line or "-----" in line:
                    if "The first frequency" in line: break
                    if "-----" in line and curr > ir_start + 10: break # End of block
                    
                parts = line.split()
                if len(parts) > 3 and ":" in parts[0]:
                    try:
                        idx = int(parts[0].replace(":", ""))
                        inten = float(parts[3])
                        if 0 <= idx < len(self.data["frequencies"]):
                            self.data["frequencies"][idx]["ir"] = inten
                    except: pass
                curr += 1

        # 3. Raman
        raman_start = -1
        for i, line in enumerate(self.lines):
            if "RAMAN SPECTRUM" in line:
                raman_start = i
                
        if raman_start != -1:
            curr = raman_start + 5
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "The first frequency" in line or "-----" in line:
                     if "The first frequency" in line: break
                     if "-----" in line and curr > raman_start + 10: break

                parts = line.split()
                if len(parts) > 2 and ":" in parts[0]:
                     try:
                        idx = int(parts[0].replace(":", ""))
                        # Mode Freq Activity ...
                        act = float(parts[2])
                        if 0 <= idx < len(self.data["frequencies"]):
                            self.data["frequencies"][idx]["raman"] = act
                     except: pass
                curr += 1

        # 4. Normal Modes
        modes_start = -1
        for i, line in enumerate(self.lines):
            if "NORMAL MODES" in line:
                modes_start = i
        
        if modes_start != -1 and self.data["atoms"]:
            n_atoms = len(self.data["atoms"])
            n_coords = n_atoms * 3
            curr = modes_start + 7
            
            mode_buffer = {} # m_idx -> [vals]
            
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line:
                    curr += 1
                    continue
                if "IR SPECTRUM" in line or "--------" in line: break
                
                try:
                    headers = [int(x) for x in line.split()]
                    start_data = curr + 1
                    for r in range(n_coords):
                        if start_data + r >= len(self.lines): break
                        dline = self.lines[start_data+r]
                        dparts = dline.split()
                        values = [float(x) for x in dparts[1:]]
                        
                        for c, m_idx in enumerate(headers):
                             if m_idx not in mode_buffer: mode_buffer[m_idx] = []
                             if c < len(values):
                                 mode_buffer[m_idx].append(values[c])
                                 
                    curr = start_data + n_coords
                except ValueError:
                    curr += 1
                    continue
            
            # Format vectors
            for m_idx, vec_flat in mode_buffer.items():
                if 0 <= m_idx < len(self.data["frequencies"]):
                    vecs = []
                    for k in range(0, len(vec_flat), 3):
                        if k+2 < len(vec_flat):
                            vecs.append((vec_flat[k], vec_flat[k+1], vec_flat[k+2]))
                    self.data["frequencies"][m_idx]["vector"] = vecs

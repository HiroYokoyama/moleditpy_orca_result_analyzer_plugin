import re
# from .logger import Logger 

class OrcaParser:
    """Parser for ORCA quantum chemistry output files"""
    def __init__(self):
        # self.logger = Logger.get_logger("OrcaParser")
        self.filename = ""
        self.raw_content = ""
        self.lines = []
        self.data = {
            "scf_traces": [], # List of SCF energy values per iteration
            "converged": False,
            "scf_energy": None,
            "atoms": [],
            "coords": [],
            "charge": 0,
            "mult": 1,
            "frequencies": [], # List of dicts: {freq, ir, raman, vector}
            "excitation_energies": [], # TDDFT
            "dipole": None,
            "dipoles": None,
            "nmr_shielding": [],
            "charges": {}, # Type -> List
            "version": None
        }

    def parse_xyz_content(self, content):
        """Parse multi-frame XYZ content."""
        lines = content.splitlines()
        steps = []
        i = 0
        n_lines = len(lines)
        
        while i < n_lines:
            line = lines[i].strip()
            if not line:
                i += 1
                continue
                
            # Number of atoms
            try:
                natoms = int(line)
            except ValueError:
                i += 1
                continue
                
            i += 1
            if i >= n_lines: break
            
            # Comment line (extract energy if possible)
            comment = lines[i].strip()
            
            # Filter out TS/CI steps if requested (User requirement: "if ts or ci present, do not use that in graph")
            # We want to filter "TS" or "CI" labels but NOT "CI-NEB" which is the method name
            import re
            # Check for TS or CI as whole words
            # Exclude strict "CI-NEB" matches from being flagged if "CI" is found
            # Simple approach: Check word boundaries. 
            # If "CI-NEB" is present, we might still have "CI" separately?
            # Let's assume labels like "Image 5 (CI)" or "TS Structure"
            
            is_excluded = False
            upper_comment = comment.upper()
            
            # Check TS
            if re.search(r"\bTS\b", upper_comment):
                is_excluded = True
                
            # Check CI, but carefully
            # Only exclude if it's NOT "CI-NEB" or similar method string
            if re.search(r"\bCI\b", upper_comment):
                 # It contains CI. Check if it is part of CI-NEB
                 if "CI-NEB" not in upper_comment:
                     is_excluded = True
            
            if is_excluded:
                # Skip this step (atoms + coords)
                # We need to advance i by natoms
                i += 1 + natoms
                continue

            energy = 0.0
            # Try to find energy in comment (e.g. "Energy: -123.456" or just "-123.456")
            # ORCA NEB format might have standard comment
            # Look for floating point number
            # Common formats: "Energy: -123.4", "Step 1 -123.4", "-123.4"
            # Let's try to find a float that looks like an energy
            floats = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+\.?", comment)
            if floats:
                try: energy = float(floats[-1]) # Take the last one? often energy is at end
                except: pass
            
            i += 1
            
            atoms = []
            coords = []
            
            for _ in range(natoms):
                if i >= n_lines: break
                parts = lines[i].split()
                if len(parts) >= 4:
                    atoms.append(parts[0])
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                i += 1
                
            steps.append({
                'type': 'neb_step',
                'energy': energy,
                'atoms': atoms,
                'coords': coords
            })
            
        return steps


    def load_from_memory(self, content, filename=""):
        self.filename = filename
        self.raw_content = content
        self.lines = content.splitlines()
        self.parse_all()

    def parse_all(self):
        self.parse_basic()
        self.parse_gradients() # Must come before trajectory to link gradients!
        self.parse_trajectory()
        self.parse_frequencies()
        self.parse_thermal()
        self.parse_mo_coeffs()
        self.parse_orbital_energies()
        self.parse_charges()
        self.parse_dipole()
        self.parse_tddft()
        self.parse_nmr()
        self.parse_basis_set()
        self.parse_scf_trace()
        
    def parse_basic(self):
        """Parse basic info: SCF Energy, Convergence, Geometry."""
        for i, line in enumerate(self.lines):
            if "Program Version" in line:
                try:
                    self.data["version"] = line.split("Version")[-1].strip().split()[0]
                except: pass

            line = line.strip()
            uu = line.upper()
            if "FINAL SINGLE POINT ENERGY" in uu:
                try:
                    self.data["scf_energy"] = float(line.split()[-1])
                except: pass
            if "TOTAL CHARGE" in uu:
                # Could be "Total Charge 0" or "Total Charge ... 0"
                try:
                    parts = line.split()
                    val = int(parts[-1])
                    self.data["charge"] = val
                except: pass
                
            if "MULTIPLICITY" in uu:
                try:
                    parts = line.split()
                    val = int(parts[-1])
                    self.data["mult"] = val
                except: pass
            if "SCF CONVERGED" in uu or "OPTIMIZATION CONVERGED" in uu or "HURRAY" in uu:
                self.data["converged"] = True
                
            if "CURRENT TRAJECTORY WILL BE WRITTEN TO" in uu:
                # Robust regex extraction
                match = re.search(r"Current trajectory will be written to\s*\.+\s*(.+)", line, re.IGNORECASE)
                if match:
                     self.data["neb_trj_file"] = match.group(1).strip()
                else:
                     # Fallback
                     try: self.data["neb_trj_file"] = line.split()[-1].strip()
                     except: pass
                
            if "CARTESIAN COORDINATES" in uu and "A.U." not in uu: # Prefer Angstrom
                # Read geometry
                self.data["atoms"] = []
                self.data["coords"] = []
                curr = i + 2
                while curr < len(self.lines):
                    l_geo = self.lines[curr].strip()
                    if not l_geo or "---" in l_geo: break
                    parts = l_geo.split()
                    if len(parts) >= 4:
                        self.data["atoms"].append(parts[0])
                        self.data["coords"].append([float(parts[1]), float(parts[2]), float(parts[3])])
                    curr += 1

    def parse_trajectory(self):
        """Parse Optimization, Scan, and NEB Trajectories."""
        self.data["scan_steps"] = []
        
        # Look for "RELAXED SURFACE SCAN STEP" or "GEOMETRY OPTIMIZATION CYCLE" or "NEB"
        # And capture Energy + Geometry
        # Actually usually Step header -> Energy -> ... -> Coordinates
        
        current_step = {}
        in_step = False
        current_scan_step = None
        
        # Helper to find coords after a header
        def read_coords_from(idx):
            atoms = []
            coords = []
            # curr = idx + 2 # Skip header and rule
            # ORCA output for coords in opt steps usually:
            # "CARTESIAN COORDINATES (ANGSTROEM)"
            # search forward for coordinates
            limit = 1000 # Search limit
            found_coords = False
            for k in range(limit):
                if idx + k >= len(self.lines): break
                line = self.lines[idx+k].strip()
                if "CARTESIAN COORDINATES" in line.upper():
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
            uu_line = line.upper()
            
            # --- NEB Parsing ---
            if "PATH SUMMARY" in uu_line and "----" in self.lines[i+1]:
                 # Found the summary table
                 # Skip header lines
                 # Line i: ---------------------------
                 # Line i+1:          PATH SUMMARY
                 # Line i+2: ---------------------------
                 # Line i+3: All forces in Eh/Bohr.
                 # Line i+4: Image Dist.(Ang.)    E(Eh)   dE(kcal/mol)  max(|Fp|)  RMS(Fp)
                 
                 curr = i + 1
                 header_found = False
                 while curr < len(self.lines) and curr < i + 10:
                     if "Image" in self.lines[curr] and "E(Eh)" in self.lines[curr]:
                         header_found = True
                         curr += 1
                         break
                     curr += 1
                 
                 if header_found:
                     while curr < len(self.lines):
                         l_row = self.lines[curr].strip()
                         if not l_row: 
                             break
                             
                         parts = l_row.split()
                         if len(parts) >= 3 and parts[0].isdigit():
                             try:
                                 img_idx = int(parts[0])
                                 dist = float(parts[1])
                                 en = float(parts[2])
                                 
                                 # We have no geometry, but user said "NO STRUCTURE IS OK"
                                 # We provide empty atoms/coords
                                 self.data["scan_steps"].append({
                                    'type': 'neb_image',
                                    'scan_step_id': img_idx,
                                    'step': len(self.data["scan_steps"]) + 1,
                                    'id': img_idx,
                                    'dist': dist,
                                    'energy': en,
                                    'atoms': [],
                                    'coords': []
                                 })
                             except: pass
                         curr += 1


            # Scan Step Header
            if "RELAXED SURFACE SCAN STEP" in uu_line:
                step_idx = 0
                match = re.search(r"STEP\s+(\d+)", uu_line)
                if match:
                    step_idx = int(match.group(1))
                current_scan_step = step_idx
                
                # Find next step to bound search
                next_marker = len(self.lines)
                for m in range(i + 1, len(self.lines)):
                    if "RELAXED SURFACE SCAN STEP" in self.lines[m].upper():
                        next_marker = m
                        break
                
                en = 0.0
                conv_info = {}
                for k in range(i, next_marker):
                    uu = self.lines[k].strip().upper()
                    if "FINAL SINGLE POINT ENERGY" in uu:
                         try: en = float(self.lines[k].split()[-1])
                         except: pass
                    elif "TOTAL ENERGY" in uu and ":" in uu and "EH" in uu:
                         # For ORCA 6: Total Energy       :        -79.79102291629319 Eh
                         try:
                             parts = self.lines[k].split(":")
                             en = float(parts[1].split()[0])
                         except: pass
                    elif "CURRENT ENERGY" in uu and "...." in uu:
                         # For ORCA relaxation blocks: Current Energy                          ....   -79.800115921 Eh
                         try:
                             parts = self.lines[k].split("....")
                             en = float(parts[1].split()[0])
                         except: pass
                    elif "GEOMETRY CONVERGENCE" in uu or "CONVERGENCE CRITERIA" in uu:
                        c_idx = k + 1
                        found_any = False
                        while c_idx < next_marker and c_idx < k + 30:
                            cl = self.lines[c_idx].strip()
                            # Only break on rule if we've already found some data lines
                            # ORCA 6.1.1 has intermediate separators, so we shouldn't break immediately.
                            if "---" in cl:
                                c_idx += 1
                                continue
                            
                            p = cl.split()
                            if len(p) >= 4:
                                s = p[-1]
                                t = p[-2]
                                v = p[-3]
                                n = " ".join(p[:-3]).strip().lower()
                                
                                # Check if it's a standard Yes/No criterion
                                if s.upper() in ["YES", "NO"]:
                                    if n and n != "item": 
                                        conv_info[n] = {
                                            'value': v,
                                            'tolerance': t,
                                            'converged': s
                                        }
                                        found_any = True
                                elif "max(" in cl.lower():
                                    # Parse Max(...) stats
                                    # e.g. Max(Bonds) 0.123 Max(Angles) 0.0
                                    matches = re.findall(r"(Max\([^)]+\))\s+([-\d\.]+)", cl, re.IGNORECASE)
                                    for label, val in matches:
                                        conv_info[label] = {
                                            'value': val, 
                                            'tolerance': '', 
                                            'converged': 'INFO'
                                        }
                                        found_any = True
                            c_idx += 1

                # Find matching gradients for this step
                step_grads = []
                candidates = []
                for g_block in self.data.get("all_gradients", []):
                    if g_block['line'] >= i and g_block['line'] < next_marker:
                        candidates.append(g_block['grads'])
                
                if candidates:
                    step_grads = candidates[-1]
                # Fallback: if not found between markers, maybe it's slightly before the marker? 
                # Or just use the one closest to the coordinate block.

                atoms, coords, found = read_coords_from(i)
                if found:
                    self.data["scan_steps"].append({
                        'type': 'scan_step',
                        'scan_step_id': current_scan_step,
                        'step': step_idx,
                        'energy': en,
                        'atoms': atoms,
                        'coords': coords,
                        'convergence': conv_info,
                        'gradients': step_grads
                    })
                    
            elif "OPTIMIZATION CYCLE" in uu_line:
                cycle_idx = 0
                match = re.search(r"CYCLE\s+(\d+)", uu_line)
                if match:
                    cycle_idx = int(match.group(1))
                
                # Find next cycle to bound search
                next_marker = len(self.lines)
                for m in range(i + 1, len(self.lines)):
                    uu_m = self.lines[m].upper()
                    # Termination markers
                    if "OPTIMIZATION CYCLE" in uu_m:
                        next_marker = m
                        break
                    if "OPTIMIZATION HAS CONVERGED" in uu_m or "OPTIMIZATION HAS RUN OUT OF CYCLES" in uu_m:
                        next_marker = m
                        break
                    if "ORCA TERMINATED NORMALLY" in uu_m:
                        next_marker = m
                        break

                en = 0.0
                conv_info = {}
                
                for k in range(i, next_marker):
                    uu = self.lines[k].strip().upper()
                    if "FINAL SINGLE POINT ENERGY" in uu:
                         try: en = float(self.lines[k].split()[-1])
                         except: pass
                    elif "TOTAL ENERGY" in uu and ":" in uu and "EH" in uu:
                         try:
                             parts = self.lines[k].split(":")
                             en = float(parts[1].split()[0])
                         except: pass
                    elif "CURRENT ENERGY" in uu and "...." in uu:
                         try:
                             parts = self.lines[k].split("....")
                             en = float(parts[1].split()[0])
                         except: pass
                    elif "GEOMETRY CONVERGENCE" in uu or "CONVERGENCE CRITERIA" in uu:
                        c_idx = k + 1
                        found_any = False
                        while c_idx < next_marker and c_idx < k + 30:
                            cl = self.lines[c_idx].strip()
                            if not cl:
                                c_idx += 1
                                continue
                            if "---" in cl:
                                c_idx += 1
                                continue
                            
                            p = cl.split()
                            if len(p) >= 4:
                                s = p[-1]
                                t = p[-2]
                                v = p[-3]
                                n = " ".join(p[:-3]).strip().lower()
                                
                                # Check if it's a standard Yes/No criterion
                                if s.upper() in ["YES", "NO"]:
                                    if n and n != "item":
                                        conv_info[n] = {
                                            'value': v,
                                            'tolerance': t,
                                            'converged': s
                                        }
                                        found_any = True
                                elif "max(" in cl.lower():
                                    # Parse Max(...) stats
                                    matches = re.findall(r"(Max\([^)]+\))\s+([-\d\.]+)", cl, re.IGNORECASE)
                                    for label, val in matches:
                                        conv_info[label] = {
                                            'value': val, 
                                            'tolerance': '', 
                                            'converged': 'INFO'
                                        }
                                        found_any = True
                            c_idx += 1
                
                # Find matching gradients for this cycle
                # Gradient block for cycle N is usually printed AFTER the convergence checks of cycle N
                step_grads = []
                # Strategy: 
                # 1. Take the LAST gradient block that appeared before next_marker
                # 2. But it MUST be at or after the current cycle index (i)
                candidates = []
                for g_block in self.data.get("all_gradients", []):
                    if g_block['line'] >= i and g_block['line'] < next_marker:
                        candidates.append(g_block['grads'])
                
                if candidates:
                    step_grads = candidates[-1] # Usually only one, but take the last if multiple
                
                # Special case: if we are at cycle N, and the gradient was printed just BEFORE the header?
                # This doesn't usually happen in ORCA, but for robustness we could check the previous few lines.
                
                # If we still don't have gradients, look slightly further back? 
                # sometimes printed just before? No, usually after.

                atoms, coords, found = read_coords_from(i)
                if found:
                     self.data["scan_steps"].append({
                        'type': 'opt_cycle',
                        'scan_step_id': current_scan_step,
                        'step': cycle_idx,
                        'energy': en,
                        'atoms': atoms,
                        'coords': coords,
                        'convergence': conv_info,
                        'gradients': step_grads
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
            uu = line.upper()
            if "MOLECULAR ORBITALS" in uu and i+1 < len(self.lines) and "---" in self.lines[i+1]:
                start_indices.append((i, "restricted"))
            elif "SPIN UP ORBITALS" in uu and i+1 < len(self.lines) and "---" in self.lines[i+1]:
                start_indices.append((i, "alpha"))
            elif "SPIN DOWN ORBITALS" in uu and i+1 < len(self.lines) and "---" in self.lines[i+1]:
                start_indices.append((i, "beta"))
                
        if not start_indices: return
        
        # Only keep the LAST occurrence for each spin type
        final_indices = {}
        for idx, spin in start_indices:
            final_indices[spin] = idx
            
        filtered_indices = []
        for spin, idx in final_indices.items():
            filtered_indices.append((idx, spin))
        filtered_indices.sort()
        
        for start_idx, spin in filtered_indices:
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
                        key = f"{idx}_{spin}"
                        if key not in self.data["mo_coeffs"]:
                             self.data["mo_coeffs"][key] = {'coeffs': [], 'spin': spin, 'energy': 0.0, 'occ': 0.0, 'id': idx}
                    
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
                                     key = f"{idx}_{spin}"
                                     if key in self.data["mo_coeffs"]:
                                         self.data["mo_coeffs"][key]['energy'] = energies[k]
                                         self.data["mo_coeffs"][key]['occ'] = occs[k]
                                 have_energy_occ = True
                                 curr += 2 # Skip these 2 lines
                             except: pass
                    except: pass
                             
                    curr += 1
                    continue
                
                # Coefficient line: "0   C  1s   0.000 ..." or "0C  1s ..."
                if len(parts) >= 2:
                     try:
                         atom_idx = -1
                         sym = ""
                         orb = ""
                         val_strs = []
                         
                         # Check for merged format "0C"
                         match_merged = re.match(r"^(\d+)([A-Za-z]+)$", parts[0])
                         if match_merged:
                             atom_idx = int(match_merged.group(1))
                             sym = match_merged.group(2)
                             orb = parts[1]
                             val_strs = parts[2:]
                         elif len(parts) >= 3 and parts[0].isdigit():
                             atom_idx = int(parts[0])
                             sym = parts[1]
                             orb = parts[2]
                             val_strs = parts[3:]
                         else:
                             curr += 1
                             continue

                         # Check if remaining parts match number of current MOs
                         if len(val_strs) == len(current_mos):
                              for k, v_str in enumerate(val_strs):
                                  mo_idx = current_mos[k]
                                  key = f"{mo_idx}_{spin}"
                                  try:
                                      val = float(v_str)
                                      if key in self.data["mo_coeffs"]:
                                          self.data["mo_coeffs"][key]['coeffs'].append({
                                              "atom_idx": atom_idx,
                                              "sym": sym,
                                              "orb": orb,
                                              "coeff": val
                                          })
                                  except: pass
                     except: pass
                curr += 1

    def parse_scan(self):
        """Alias for parse_trajectory."""
        self.parse_trajectory()

    def parse_gradient(self):
        """Alias for parse_gradients."""
        self.parse_gradients()

    def parse_gradients(self):
        """Parse all Cartesian Gradient blocks found in the file."""
        self.data["gradients"] = [] # The last one (default)
        self.data["all_gradients"] = [] # List of {line: int, grads: []}
        
        gradient_starts = []
        for i, line in enumerate(self.lines):
            stripped = line.strip().upper()
            if "CARTESIAN GRADIENT" in stripped and "NORM" not in stripped:
                gradient_starts.append(i)
        
        if not gradient_starts: return
        
        for start_idx in gradient_starts:
            block_grads = []
            curr = start_idx + 1
            # Skip header separators or empty lines until data matches format
            found_data = False
            while curr < len(self.lines) and curr < start_idx + 15:
                line = self.lines[curr].strip()
                parts = line.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    found_data = True
                    break 
                curr += 1
            
            if not found_data: continue

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "-------" in line or "Difference to" in line: break
                if not line: break 
                
                parts = line.split()
                # Format: 0 C : X Y Z (6 parts) or 1 C X Y Z (5 parts)
                if len(parts) >= 5:
                    try:
                        if parts[2] == ":":
                            if len(parts) >= 6:
                                idx_raw = int(parts[0])
                                sym = parts[1]
                                vx, vy, vz = float(parts[3]), float(parts[4]), float(parts[5])
                                # Usually 1-indexed
                                idx = idx_raw - 1 if idx_raw > 0 else 0
                            else: 
                                curr += 1
                                continue
                        else:
                            idx_raw = int(parts[0])
                            sym = parts[1]
                            vx, vy, vz = float(parts[2]), float(parts[3]), float(parts[4])
                            idx = idx_raw - 1 if idx_raw > 0 else 0
                            
                        block_grads.append({
                            'atom_idx': idx, 
                            'atom_sym': sym, 
                            'vector': [vx, vy, vz]
                        })
                    except: pass
                curr += 1
            
            if block_grads:
                self.data["all_gradients"].append({
                    'line': start_idx,
                    'grads': block_grads
                })
        
        if self.data["all_gradients"]:
            # Set the last one as the default "gradients"
            self.data["gradients"] = self.data["all_gradients"][-1]['grads']

    def parse_dipole(self):
        # Look for "Total Dipole Moment"
        self.data["dipoles"] = None
        self.data["dipole"] = None
        
        candidates = []
        for i, line in enumerate(self.lines):
            uu = line.upper()
            if "TOTAL DIPOLE MOMENT" in uu and ":" in line:
                 candidates.append(i)
                 
        if not candidates: return
        
        # Take the last occurrence
        idx = candidates[-1]
        line = self.lines[idx]
        parts = line.split(":")
        if len(parts) > 1:
            try:
                vec_str = parts[1].strip().split()
                if len(vec_str) >= 3:
                     x, y, z = float(vec_str[0]), float(vec_str[1]), float(vec_str[2])
                     
                     mag = 0.0
                     if idx + 1 < len(self.lines):
                         line2 = self.lines[idx+1]
                         if "Magnitude" in line2 and ":" in line2:
                             mag = float(line2.split(":")[1].strip())
                         else:
                             import math
                             mag = math.sqrt(x*x + y*y + z*z)
                             
                     self.data["dipoles"] = {
                         "vector": (x, y, z),
                         "magnitude": mag
                     }
                     self.data["dipole"] = self.data["dipoles"]
            except: pass


        
    def parse_charges(self):
        self.data["charges"] = {} # type -> list of {atom_idx, atom_sym, charge}
        
        # Section Markers
        mulliken_start = -1
        loewdin_start = -1
        hirshfeld_start = -1
        mayer_start = -1
        nbo_start = -1
        chelpg_start = -1
        mk_start = -1
        mbis_start = -1
        resp_start = -1
        fmo_start = -1
        
        for i, line in enumerate(self.lines):
            uu = line.upper()
            if "MULLIKEN ATOMIC CHARGES" in uu: mulliken_start = i
            elif "LOEWDIN ATOMIC CHARGES" in uu: loewdin_start = i
            elif "HIRSHFELD ANALYSIS" in uu: hirshfeld_start = i
            elif "MAYER POPULATION ANALYSIS" in uu: mayer_start = i
            elif "NATURAL POPULATIONS" in uu: nbo_start = i
            elif "CHELPG ATOMIC CHARGES" in uu: chelpg_start = i
            elif "MERZ-KOLLMAN ATOMIC CHARGES" in uu or "MK ATOMIC CHARGES" in uu: mk_start = i
            elif "MBIS ANALYSIS" in uu: mbis_start = i
            elif "RESP ATOMIC CHARGES" in uu: resp_start = i
            elif "FRONTIER MOLECULAR ORBITAL POPULATION ANALYSIS" in uu: fmo_start = i
            
        def parse_standard_block(start_idx, header_lines=2, hirshfeld=False, mbis=False):
            res = []
            if start_idx == -1: return res
            curr = start_idx + header_lines
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line or "---" in line or "Sum of" in line:
                    if res: break
                    curr += 1
                    continue
                
                parts = line.split()
                if len(parts) >= 3:
                     try:
                         # Handle merged "0C" case
                         first_part = parts[0]
                         match = re.match(r"^(\d+)([a-zA-Z]+)$", first_part)
                         if match:
                             idx_str = match.group(1)
                             sym = match.group(2)
                         else:
                             idx_str = first_part.strip(":")
                             # If "0 C", sym is parts[1]
                             sym = parts[1].strip(":")
                         
                         if hirshfeld or mbis:
                              val = float(parts[2]) # 0 C Charge Spin
                         else:
                              # 0 C : -1.23 or 0 C -1.23
                              val = float(parts[3]) if parts[2] == ":" else float(parts[2])
                         
                         atom_data = {
                             "atom_idx": int(idx_str),
                             "atom_sym": sym,
                             "charge": val
                         }
                         
                         # Capture spin if available (Hirshfeld/MBIS)
                         if hirshfeld:
                             if len(parts) >= 4:
                                 atom_data["spin"] = float(parts[3])
                         elif mbis:
                             if len(parts) >= 4:
                                 # ATOM CHARGE POPULATION [SPIN]
                                 atom_data["population"] = float(parts[3])
                             if len(parts) >= 5:
                                 atom_data["spin"] = float(parts[4])
                             
                         res.append(atom_data)
                     except: pass
                curr += 1
            return res

        self.data["charges"]["Mulliken"] = parse_standard_block(mulliken_start)
        self.data["charges"]["Loewdin"] = parse_standard_block(loewdin_start)
        self.data["charges"]["Hirshfeld"] = parse_standard_block(hirshfeld_start, hirshfeld=True)
        self.data["charges"]["CHELPG"] = parse_standard_block(chelpg_start)
        self.data["charges"]["MK"] = parse_standard_block(mk_start)
        self.data["charges"]["MBIS"] = parse_standard_block(mbis_start, header_lines=3, mbis=True)
        self.data["charges"]["RESP"] = parse_standard_block(resp_start)
        
        # Clean up empty entries
        self.data["charges"] = {k: v for k, v in self.data["charges"].items() if v}
        
        # Mayer Parsing (QA is Mulliken charge)
        if mayer_start != -1:
            mayer_res = []
            curr = mayer_start + 1
            while curr < len(self.lines) and curr < mayer_start + 15:
                if "ATOM" in self.lines[curr] and "QA" in self.lines[curr]:
                    curr += 1
                    break
                curr += 1
            
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line or "---" in line or "Mayer bond" in line: break
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        idx = int(parts[0])
                        sym = parts[1]
                        # NA, ZA, QA, VA, BVA, FA
                        # parts[2] = NA
                        # parts[3] = ZA (Atomic Number)
                        # parts[4] = QA (Charge)
                        qa = float(parts[4])
                        
                        # Extra Mayer Indices
                        extra = {}
                        if len(parts) >= 8:
                            extra["valency"] = float(parts[5])       # VA
                            extra["bonded_valency"] = float(parts[6]) # BVA
                            extra["free_valency"] = float(parts[7])   # FA
                            
                        mayer_res.append({"atom_idx": idx, "atom_sym": sym, "charge": qa, **extra})
                    except: pass
                curr += 1
            if mayer_res: 
                self.data["charges"]["Mayer"] = mayer_res
                if not self.data["charges"].get("Mulliken"):
                    self.data["charges"]["Mulliken"] = mayer_res

        # NBO Parsing
        if nbo_start != -1:
            nbo_charges = []
            
            # 1. Try to find Detailed Summary Table first
            # "Summary of Natural Population Analysis"
            summary_start = -1
            # Search a bit deeper than nbo_start, as "NATURAL POPULATIONS" might be just one header
            # But usually the summary is a distinct block.
            # Let's search the whole file? No, usually near the end or after NBO analysis
            # Optimization: Look forward from nbo_start
            
            for i in range(nbo_start, min(nbo_start + 2000, len(self.lines))):
                 if "Summary of Natural Population Analysis" in self.lines[i]:
                     summary_start = i
                     break
            
            if summary_start != -1:
                # Parse the detailed table
                curr = summary_start + 1
                while curr < len(self.lines):
                    line = self.lines[curr].strip()
                    if "Atom No" in line and "Charge" in line and "Core" in line:
                        curr += 1 # Skip header line
                        if curr < len(self.lines) and "----" in self.lines[curr]: 
                            curr += 1 # Skip separator
                        break
                    curr += 1
                
                while curr < len(self.lines):
                    line = self.lines[curr].strip()
                    if "====" in line or "Total" in line: break
                    if not line:
                        curr += 1
                        continue
                    
                    parts = line.split()
                    # C  1   -0.49495      1.99999     4.48487    0.01009     6.49495
                    if len(parts) >= 7:
                        try:
                            sym = parts[0]
                            idx = int(parts[1])
                            chg = float(parts[2])
                            core = float(parts[3])
                            val = float(parts[4])
                            ryd = float(parts[5])
                            tot = float(parts[6])
                            
                            nbo_charges.append({
                                "atom_idx": idx - 1,
                                "atom_sym": sym,
                                "charge": chg,
                                "core": core,
                                "valence": val,
                                "rydberg": ryd,
                                "total": tot
                            })
                        except: pass
                    curr += 1
            
            # 2. Fallback if no summary found (or parsing failed), try simple block near nbo_start
            if not nbo_charges:
                curr = nbo_start + 1
                while curr < len(self.lines):
                    line = self.lines[curr].strip()
                    if "---" in line:
                        curr += 1
                        continue
                    if "================" in line or "Natural Electron Configuration" in line:
                        if nbo_charges: break
                        curr += 1
                        continue
                    if not line:
                        if nbo_charges: break # Stop on empty line if we have data
                        curr += 1
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 3:
                         try:
                             # Simple format: Atom No Charge ...
                             # or format: C 1 -0.123
                             sym = parts[0]
                             idx = int(parts[1])
                             chg = float(parts[2])
                             nbo_charges.append({
                                 "atom_idx": idx - 1, 
                                 "atom_sym": sym,
                                 "charge": chg
                             })
                         except: pass
                    curr += 1

            if nbo_charges:
                self.data["charges"]["NBO"] = nbo_charges

        # FMO Parsing
        if fmo_start != -1:
             fmo_data = []
             curr = fmo_start + 1
             
             table_start = False
             # Look for table header
             while curr < len(self.lines) and curr < fmo_start + 40:
                 if "--------" in self.lines[curr]:
                      # Check previous lines for "Atom" "HOMO" etc.
                      # Usually two or three lines of header, then barrier
                      # Or "Atom   Q(Mulliken) ..."
                      if curr > 0 and "Atom" in self.lines[curr-1]:
                          table_start = True
                          curr += 1
                          break
                      if curr > 1 and "HOMO" in self.lines[curr-2]:
                          table_start = True
                          curr += 1
                          break
                 curr += 1
            
             if table_start:
                 while curr < len(self.lines):
                     line = self.lines[curr].strip()
                     if not line or "--------" in line:
                         if fmo_data: break
                         curr += 1
                         continue
                     
                     parts = line.split()
                     # 0-C 0.937 0.906 0.804 0.755
                     if len(parts) >= 5:
                         try:
                             atom_lbl = parts[0] # 0-C
                             if "-" in atom_lbl:
                                 # Split only on first hyphen to handle negative indices or strange names?
                                 # Usually ORCA format is 0-C, 1-H etc.
                                 p_lbl = atom_lbl.split("-")
                                 idx_str = p_lbl[0] 
                                 sym = p_lbl[1]
                                 idx = int(idx_str)
                             else:
                                 # Fallback if just C or just 0
                                 idx = len(fmo_data)
                                 sym = atom_lbl
                             
                             homo_m = float(parts[1])
                             homo_l = float(parts[2])
                             lumo_m = float(parts[3])
                             lumo_l = float(parts[4])
                             
                             fmo_data.append({
                                 "atom_idx": idx,
                                 "atom_sym": sym,
                                 # Use Mulliken HOMO as primary visual if asked for "charge"
                                 "charge": homo_m, 
                                 "homo_mulliken": homo_m,
                                 "homo_loewdin": homo_l,
                                 "lumo_mulliken": lumo_m,
                                 "lumo_loewdin": lumo_l
                             })
                         except: pass
                     curr += 1
                 
                 if fmo_data:
                     self.data["charges"]["FMO"] = fmo_data

        
    def parse_nmr(self):
        self.data["nmr_shielding"] = [] # List of {atom_idx, atom_sym, shielding, shift=None}
        
        # Look for "CHEMICAL SHIELDING SUMMARY (ppm)" or individual nucleus blocks
        # ORCA output typically has a summary at the end of properties
        
        summary_start = -1
        # Try to find the summary table first
        for i, line in enumerate(self.lines):
            if "CHEMICAL SHIELDING SUMMARY (PPM)" in line.upper():
                summary_start = i
                
        if summary_start != -1:
            curr = summary_start + 1
            # Skip until we hit the header "N Nucleus Shielding"
            while curr < len(self.lines) and curr < summary_start + 10:
                if "N" in self.lines[curr] and "Shielding" in self.lines[curr]:
                    curr += 1
                    break
                curr += 1

            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line: break
                if "--------------" in line:
                    if len(self.data["nmr_shielding"]) > 0: break
                    curr += 1
                    continue
                
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
            uu = line.upper()
            if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in uu:
                start_idx = i
                found_block = True
                break
            # Fallback for simpler header
            if "ABSORPTION SPECTRUM" in uu and "ELECTRIC DIPOLE" in uu and not found_block:
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
        

    def parse_thermal(self):
        self.data["thermal"] = {}
        # ORCA Thermochem block
        # Look for "THERMOCHEMISTRY AT 298.15 K"
    
        start_line = -1
        # Scan for the start of thermochemistry section
        for i, line in enumerate(self.lines):
            line_upper = line.upper()
            if "THERMOCHEMISTRY AT" in line_upper:
                start_line = i
                
        if start_line != -1:
            curr = start_line
            
            # Mapping of ORCA labels to JSON keys
            # Using uppercase for case-insensitive matching
            thermal_keys = {
                "ELECTRONIC ENERGY": "electronic_energy",
                "ZERO POINT ENERGY": "zpe",
                "THERMAL VIBRATIONAL CORRECTION": "corr_vib",
                "THERMAL ROTATIONAL CORRECTION": "corr_rot",
                "THERMAL TRANSLATIONAL CORRECTION": "corr_trans",
                "TOTAL THERMAL ENERGY": "thermal_energy",
                "TOTAL THERMAL CORRECTION": "corr_thermal_total",
                "NON-THERMAL (ZPE) CORRECTION": "corr_zpe",
                "TOTAL CORRECTION": "corr_total",
                "TOTAL ENTHALPY": "enthalpy",
                "THERMAL ENTHALPY CORRECTION": "thermal_enthalpy_corr", # Just RT
                "ELECTRONIC ENTROPY": "s_el",
                "VIBRATIONAL ENTROPY": "s_vib",
                "ROTATIONAL ENTROPY": "s_rot",
                "TRANSLATIONAL ENTROPY": "s_trans",
                "FINAL ENTROPY TERM": "entropy", 
                "FINAL GIBBS FREE ENERGY": "gibbs",
                "G-E(EL)": "gibbs_corr"
            }
            
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                line_upper = line.upper()
                
                if "TIMINGS FOR INDIVIDUAL MODULES" in line_upper: break
                
                # Temperature check
                if "TEMPERATURE" in line_upper and "K" in line_upper and "..." in line_upper:
                    match = re.search(r"(\d+\.\d+)\s*K", line)
                    if match:
                         self.data["thermal"]["temperature"] = float(match.group(1))

                # Iterate through expected keys
                for key_upper, val_key in thermal_keys.items():
                    if key_upper in line_upper:
                        # Extract the numeric value, prioritizing heartrees (Eh)
                        if "Eh" in line:
                            # Split by "Eh" and take the part before it
                            pre_eh = line.split("Eh")[0]
                            # Find the last float-like number in the pre-Eh part
                            val_matches = re.findall(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?", pre_eh)
                            if val_matches:
                                try:
                                    val = float(val_matches[-1])
                                    self.data["thermal"][val_key] = val
                                except ValueError: pass
                        else:
                            # Fallback: take the last numeric match
                            matches = re.findall(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?", line)
                            if matches:
                                try:
                                    val = float(matches[-1])
                                    self.data["thermal"][val_key] = val
                                except ValueError: pass
                curr += 1
            
            # Post-processing: Calculate H correction (H - E_el) more accurately
            # ORCA's "Thermal Enthalpy correction" is often just RT.
            # We want the total correction including ZPE and Thermal effects.
            t_data = self.data["thermal"]
            if "enthalpy" in t_data and "electronic_energy" in t_data:
                t_data["enthalpy_corr"] = t_data["enthalpy"] - t_data["electronic_energy"]
            elif "corr_total" in t_data and "thermal_enthalpy_corr" in t_data:
                # Fallback: H_corr = Total_corr (U-E_el) + RT_corr
                t_data["enthalpy_corr"] = t_data["corr_total"] + t_data["thermal_enthalpy_corr"]
            
            # Similarly for Gibbs if not directly found
            if "gibbs" in t_data and "electronic_energy" in t_data:
                 t_data["gibbs_corr"] = t_data["gibbs"] - t_data["electronic_energy"]
            
            # Count imaginary frequencies
            freqs = self.data.get("frequencies", [])
            imaginary_count = sum(1 for f in freqs if f.get("freq", 0) < 0)
            t_data["imaginary_freq_count"] = imaginary_count
            
            return 

            
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



    def parse_frequencies(self):
        self.data["frequencies"] = []
        
        # 1. Frequencies
        freq_start = -1
        for i, line in enumerate(self.lines):
            if "VIBRATIONAL FREQUENCIES" in line.upper():
                freq_start = i
                # Don't break immediately, could be multiple? usually last one matters or first?
                # In optimization + freq, it's at end.
        
        if freq_start != -1:
            curr = freq_start + 1
            # Skip until we find first data line (Index: Value)
            while curr < len(self.lines) and curr < freq_start + 10:
                line = self.lines[curr].strip()
                if ":" in line and ("cm**-1" in line or "cm-1" in line):
                    break # Found data match
                curr += 1
                
            # Now parse data
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "NORMAL MODES" in line or "-------" in line:
                    # If we already have freqs, a dashed line means end
                    if len(self.data["frequencies"]) > 0: break
                
                if not line:
                    curr += 1
                    continue
                
                # Format:   0:         0.00 cm**-1
                if ":" in line and ("cm**-1" in line or "cm-1" in line):
                    parts = line.split()
                    # Find value after colon?
                    # usually parts: "0:", "0.00", "cm**-1"
                    try:
                        # Value is typically index 1 if index 0 ends with colon
                        val_str = parts[1]
                        val = float(val_str)
                        self.data["frequencies"].append({"freq": val, "ir": 0.0, "raman": 0.0, "vector": []})
                    except: pass
                elif len(self.data["frequencies"]) > 0 and ":" not in line:
                    # Maybe end of block
                    break
                    
                curr += 1

        # 2. IR Spectrum
        ir_start = -1
        for i, line in enumerate(self.lines):
            if "IR SPECTRUM" in line.upper():
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
            if "RAMAN SPECTRUM" in line.upper():
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
            uu = line.upper()
            if "NORMAL MODES" in uu:
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
            
            # Process collected mode buffer AFTER parsing all mode blocks
            for m_idx, vec_flat in mode_buffer.items():
                if 0 <= m_idx < len(self.data["frequencies"]):
                    vecs = []
                    for k in range(0, len(vec_flat), 3):
                        if k+2 < len(vec_flat):
                            vecs.append((vec_flat[k], vec_flat[k+1], vec_flat[k+2]))
                    self.data["frequencies"][m_idx]["vector"] = vecs



    def parse_orbital_energies(self):
        """Parse orbital energies from ORBITAL ENERGIES section"""
        self.data["orbital_energies"] = []
        self.data["mos"] = [] # Backward compatibility
        
        # Look for "ORBITAL ENERGIES" section
        start_indices = []
        for i, line in enumerate(self.lines):
            uu = line.upper()
            if "ORBITAL ENERGIES" in uu and i+1 < len(self.lines) and "---" in self.lines[i+1]:
                start_indices.append((i, "restricted")) # Default
            elif "SPIN UP ORBITALS" in uu and i+1 < len(self.lines) and "---" in self.lines[i+1]:
                start_indices.append((i, "alpha"))
            elif "SPIN DOWN ORBITALS" in uu and i+1 < len(self.lines) and "---" in self.lines[i+1]:
                start_indices.append((i, "beta"))
        
        if not start_indices:
            return
            
        # Only keep the LAST occurrence for each spin type
        final_indices = {}
        for idx, spin in start_indices:
            final_indices[spin] = idx
            
        filtered_indices = []
        for spin, idx in final_indices.items():
            filtered_indices.append((idx, spin))
        filtered_indices.sort()
        
        for start_idx, spin in filtered_indices:
            # Skip title and separator, then find column header
            curr = start_idx + 2
            while curr < len(self.lines) and curr < start_idx + 10:
                line = self.lines[curr].strip()
                # Check for "NO", "OCC", and "Eh" or "eV"
                if "NO" in line and ("OCC" in line or "E(Eh)" in line or "E(eV)" in line):
                    curr += 1
                    break
                curr += 1
                
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if not line or "---" in line or "****" in line or "MULLIKEN" in line:
                    break
                
                parts = line.split()
                # Format: NO   OCC          E(Eh)            E(eV)
                #           0   2.0000     -10.186408      -277.1862
                if len(parts) >= 4:
                    try:
                        orbital_idx = int(parts[0])
                        occupation = float(parts[1])
                        energy_eh = float(parts[2])
                        energy_ev = float(parts[3])
                        
                        orb_data = {
                            "index": orbital_idx,
                            "id": orbital_idx, # For mo_analysis
                            "occupation": occupation,
                            "occ": occupation, # For mo_analysis
                            "energy_eh": energy_eh,
                            "energy_ev": energy_ev,
                            "energy": energy_eh, # For backward compatibility
                            "spin": spin,
                            "type": "occupied" if occupation > 0.1 else "virtual"
                        }
                        self.data["orbital_energies"].append(orb_data)
                        self.data["mos"].append(orb_data)
                    except:
                        pass
                
                curr += 1
    
    # COMMENTED OUT: Hessian file parsing no longer needed
    # Force analysis now uses only output file gradient data
    # def parse_hessian_file(self, filepath):
    #     """Parse ORCA Hessian (.hess) file to extract force constants"""
    #     try:
    #         with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
    #             content = f.read()
    #     except Exception as e:
    #         print(f"Error reading Hessian file: {e}")
    #         return None
    #     
    #     lines = content.splitlines()
    #     hessian_data = {
    #         "matrix": [],
    #         "frequencies": [],
    #         "normal_modes": [],
    #         "atoms": [],
    #         "coords": [],
    #         "gradient_vec": []
    #     }
    #     
    #     i = 0
    #     while i < len(lines):
    #         line = lines[i].strip()
    #         
    #         # Parse Hessian matrix
    #         if line == "$hessian":
    #             i += 1
    #             if i < len(lines):
    #                 n_coords = int(lines[i].strip())
    #                 i += 1
    #                 
    #                 # Read matrix in blocks
    #                 matrix = [[] for _ in range(n_coords)]
    #                 
    #                 while i < len(lines):
    #                     line = lines[i].strip()
    #                     if line.startswith("$"):
    #                         break
    #                     
    #                     parts = line.split()
    #                     if len(parts) > 1 and parts[0].isdigit():
    #                         # Check if it's a header line (all integers)
    #                         is_header = True
    #                         try:
    #                             for p in parts:
    #                                 if "." in p: # Floats usually have decimal points in ORCA
    #                                     is_header = False
    #                                     break
    #                                 # Attempt to parse as float? ORCA floats are 1.23E-01
    #                                 # Integers don't have '.' usually
    #                         except:
    #                             pass
    #                             
    #                         if is_header:
    #                             # Double check: try parsing all as ints
    #                             try:
    #                                 [int(x) for x in parts]
    #                                 # It is a header, skip it
    #                                 i += 1
    #                                 continue
    #                             except ValueError:
    #                                 # Not all ints, so it's data
    #                                 pass
    #
    #                         row_idx = int(parts[0])
    #                         values = [float(x) for x in parts[1:]]
    #                         if row_idx < len(matrix):
    #                             matrix[row_idx].extend(values)
    #                     
    #                     i += 1
    #                 
    #                 hessian_data["matrix"] = matrix
    #         
    #         # Parse frequencies
    #         elif line == "$vibrational_frequencies":
    #             i += 1
    #             if i < len(lines):
    #                 n_freq = int(lines[i].strip())
    #                 i += 1
    #                 
    #                 frequencies = []
    #                 for _ in range(n_freq):
    #                     if i >= len(lines):
    #                         break
    #                     parts = lines[i].strip().split()
    #                     if len(parts) >= 2:
    #                         frequencies.append(float(parts[1]))
    #                     i += 1
    #                 
    #                 hessian_data["frequencies"] = frequencies
    #         
    #         # Parse normal modes
    #         elif line == "$normal_modes":
    #             i += 1
    #             if i < len(lines):
    #                 dims = lines[i].strip().split() # rows cols
    #                 if len(dims) >= 2:
    #                     n_rows = int(dims[0]) # 3 * n_atoms
    #                     n_cols = int(dims[1]) # n_freq
    #                     i += 1
    #                     
    #                     # Initialize normal modes matrix (n_rows x n_cols)
    #                     # We want to store per mode, so list of n_cols vectors of length n_rows
    #                     modes = [[] for _ in range(n_cols)]
    #                     
    #                     while i < len(lines):
    #                         line = lines[i].strip()
    #                         if line.startswith("$"):
    #                             break
    #                         
    #                         # Block headers: 0 1 2 3 ...
    #                         parts = line.split()
    #                         if not parts:
    #                             i += 1
    #                             continue
    #                             
    #                         # Check if header line (likely integer indices)
    #                         try:
    #                             # If all are integers, it's a header line
    #                             col_indices = [int(x) for x in parts]
    #                             i += 1 # Move to data
    #                             
    #                             # Read data rows for this block
    #                             for r in range(n_rows):
    #                                 if i >= len(lines): break
    #                                 line = lines[i].strip()
    #                                 parts = line.split()
    #                                 # First part is row index
    #                                 if len(parts) > 1:
    #                                     # row_idx = int(parts[0])
    #                                     vals = [float(x) for x in parts[1:]]
    #                                     for k, col_idx in enumerate(col_indices):
    #                                         if col_idx < len(modes):
    #                                             modes[col_idx].append(vals[k])
    #                                 i += 1
    #                         except ValueError:
    #                             # Not a header line? Should not happen if format is strict
    #                             i += 1
    #                     
    #                     # Reformat to list of vectors (tuples) for each mode
    #                     # Each mode is list of 3N floats. Convert to list of (x,y,z)
    #                     formatted_modes = []
    #                     for m in modes:
    #                         vecs = []
    #                         for k in range(0, len(m), 3):
    #                             if k+2 < len(m):
    #                                 vecs.append((m[k], m[k+1], m[k+2]))
    #                         formatted_modes.append(vecs)
    #                         
    #                     hessian_data["normal_modes"] = formatted_modes
    #
    #         # Parse atoms
    #         elif line == "$atoms":
    #             i += 1
    #             if i < len(lines):
    #                 n_atoms = int(lines[i].strip())
    #                 i += 1
    #                 
    #                 atoms = []
    #                 coords = []
    #                 for _ in range(n_atoms):
    #                     if i >= len(lines):
    #                         break
    #                     parts = lines[i].strip().split()
    #                     if len(parts) >= 5:
    #                         atoms.append(parts[0])
    #                         # Convert Bohr to Angstrom
    #                         BOHR_TO_ANG = 0.529177249
    #                         coords.append([
    #                             float(parts[2]) * BOHR_TO_ANG,
    #                             float(parts[3]) * BOHR_TO_ANG,
    #                             float(parts[4]) * BOHR_TO_ANG
    #                         ])
    #                     i += 1
    #                 
    #                 hessian_data["atoms"] = atoms
    #                 hessian_data["coords"] = coords
    #
    #         # Parse Gradient (Forces)
    #         elif line == "$gradient":
    #             i += 1
    #             if i < len(lines):
    #                 try:
    #                     n_grad = int(lines[i].strip())
    #                     i += 1
    #                     
    #                     gradients_vec = []
    #                     for _ in range(n_grad):
    #                         if i >= len(lines): break
    #                         
    #                         parts = lines[i].strip().split()
    #                         if not parts:
    #                             i += 1
    #                             continue
    #                             
    #                         val = 0.0
    #                         # Robust parsing: Check format "index value" vs "value"
    #                         if len(parts) >= 2:
    #                             # Assume "index value" format usually
    #                             # Check if first part is int
    #                             try:
    #                                 # If parts[0] is int, parts[1] is value
    #                                 idx = int(parts[0])
    #                                 val = float(parts[1])
    #                             except ValueError:
    #                                 # Maybe format is "value value"? Unlikely for gradient (1D)
    #                                 # Or parts[0] is actually the float value?
    #                                 try:
    #                                     val = float(parts[0])
    #                                 except:
    #                                     pass
    #                         elif len(parts) == 1:
    #                             # Just value
    #                             try:
    #                                 val = float(parts[0])
    #                             except: pass
    #                             
    #                         gradients_vec.append(val)
    #                         i += 1
    #                     
    #                     hessian_data["gradient_vec"] = gradients_vec
    #                 except ValueError:
    #                      print("Error parsing number of gradients")
    #                      i += 1
    #         
    #         i += 1
    #     
    #     return hessian_data
        
        
    def parse_basis_set(self):
        """Parse Basis Set information needed for MO visualization"""
        self.data["basis_set_shells"] = []
        
        # Look for "BASIS SET IN INPUT FORMAT"
        start_idx = -1
        for i, line in enumerate(self.lines):
            if "BASIS SET IN INPUT FORMAT" in line:
                start_idx = i
                break
        
        if start_idx == -1: return
        
        curr = start_idx + 2
        
        basis_defs = {} # Sym -> List of shells
        current_sym = None
        current_shells = []
        
        while curr < len(self.lines):
            line = self.lines[curr].strip()
            
            # Stop conditions
            if "--------" in line and curr > start_idx + 10: break 
            if "AUXILIARY BASIS" in line: break
            
            # Start of Atom block: "NewGTO H"
            if line.startswith("NewGTO"):
                parts = line.split()
                if len(parts) >= 2:
                    current_sym = parts[1]
                    current_shells = []
                curr += 1
                continue
                
            # End of Atom block: "end" or "end;"
            if line.startswith("end"):
                if current_sym:
                     basis_defs[current_sym] = current_shells
                curr += 1
                continue
                
            # Shell definition header: "S   3" or "P   2"
            parts = line.split()
            if len(parts) >= 2 and parts[0].upper() in ['S', 'P', 'D', 'F', 'G']:
                sh_type = parts[0].upper()
                try:
                    n_prim = int(parts[1])
                    if n_prim > 50: # Sanity check
                        curr += 1
                        continue
                        
                    curr += 1
                    
                    exps = []
                    coeffs = []
                    
                    # Read primitives
                    for _ in range(n_prim):
                        if curr >= len(self.lines): break
                        pl = self.lines[curr].strip()
                        pp = pl.split()
                        if len(pp) >= 3:
                            exps.append(float(pp[1]))
                            coeffs.append(float(pp[2]))
                        curr += 1
                        
                    l_map = {'S':0, 'P':1, 'D':2, 'F':3, 'G':4}
                    l_val = l_map.get(sh_type, 0)
                    
                    if exps:
                        current_shells.append({
                            'l': l_val,
                            'exps': exps,
                            'coeffs': coeffs
                        })
                    
                    continue 
                except: 
                    pass
            
            curr += 1
            if curr > start_idx + 5000: break # Safety break
            
        # Expand to actual atoms
        atoms = self.data.get("atoms", [])
        coords = self.data.get("coords", [])
        
        full_shells = []
        
        if not atoms:
             # Try to recover atoms from parser data if parse_basic hasn't run or failed?
             # No, parse_basis_set relies on atoms being parsed.
             return

        for idx, (sym, coord) in enumerate(zip(atoms, coords)):
            defs = basis_defs.get(sym, [])
            for d in defs:
                full_shells.append({
                    'atom_idx': idx,
                    'origin': coord,
                    'l': d['l'],
                    'exps': d['exps'],
                    'coeffs': d['coeffs']
                })
                
        self.data["basis_set_shells"] = full_shells

    def parse_scf_trace(self):
        """Extract SCF iteration energies for each block found."""
        self.data["scf_traces"] = []
        
        # We search for blocks starting with D-I-I-S or S-O-S-C-F or SCF ITERATIONS
        current_step_label = "Initial"
        
        i = 0
        while i < len(self.lines):
            line = self.lines[i]
            uu = line.upper()
        
            if "OPTIMIZATION CYCLE" in uu:
                try: 
                    parts = line.split()
                    # Find the index of CYCLE and take the next part
                    cycle_part = parts[parts.index("CYCLE") + 1] if "CYCLE" in parts else parts[-2]
                    current_step_label = f"Cycle {cycle_part}"
                except: current_step_label = "Opt Cycle"
            elif "SCAN STEP" in uu:
                try: current_step_label = f"Scan Step {line.split()[-1]}"
                except: current_step_label = "Scan Step"
            elif "ORCA PROPERTIES" in uu or "ORCA PROPERTY" in uu:
                current_step_label = "Property/Final"
            elif "OPTIMIZATION HAS CONVERGED" in uu:
                current_step_label = "Post-Opt/Final"
                
            if "SCF ITERATIONS" in uu or "ORCA LEAN-SCF" in uu or "INCREMENTAL FOCK MATRIX" in uu:
                header_idx = -1
                for k in range(1, 15):
                    if i + k >= len(self.lines): break
                    uu_k = self.lines[i+k].upper()
                    if "ITER" in uu_k and "ENERGY" in uu_k:
                        header_idx = i + k
                        break
                
                if header_idx != -1:
                    trace = []
                    idx = header_idx + 1
                    # Check for separator line right after header
                    if idx < len(self.lines) and "---" in self.lines[idx]:
                        idx += 1
                    
                    while idx < len(self.lines):
                        l_scf = self.lines[idx].strip()
                        if not l_scf or "---" in l_scf or "SUCCESS" in l_scf or "Energy Check" in l_scf:
                            if trace: break
                            idx += 1
                            continue
                        
                        parts = l_scf.split()
                        if len(parts) >= 2:
                            try:
                                it_no = int(parts[0])
                                it_en = float(parts[1])
                                trace.append({'iter': it_no, 'energy': it_en})
                            except: pass
                        idx += 1
                    
                    if trace:
                        # Check if we should append or start new
                        # If it's a new "SCF ITERATIONS" block, we want a new entry in the traces list
                        # unless it's genuinely part of the same convergence process (rare in ORCA output stream)
                        
                        # Heuristic: if last trace in self.data["scf_traces"] has same label, 
                        # suffix the new one if they are distinct blocks.
                        same_label_count = 0
                        for t in self.data["scf_traces"]:
                            if t['step'].startswith(current_step_label):
                                same_label_count += 1
                        
                        label = current_step_label
                        if same_label_count > 0:
                            label = f"{current_step_label} ({same_label_count + 1})"
                            
                        self.data["scf_traces"].append({
                            'step': label,
                            'iterations': trace
                        })
                        i = idx
                        continue
            i += 1
        # print(f"OrcaParser: Parsed {len(self.data['scf_traces'])} SCF trace blocks.")

import re
# from .logger import Logger 

class OrcaParser:
    """Parser for ORCA quantum chemistry output files"""
    def __init__(self):
        # self.logger = Logger.get_logger("OrcaParser")
        # print("OrcaParser initialized")
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
            "charges": {} # Type -> List
        }

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
        self.parse_mo()
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
            line = line.strip()
            uu = line.upper()
            if "FINAL SINGLE POINT ENERGY" in uu:
                try:
                    self.data["scf_energy"] = float(line.split()[-1])
                except: pass
            if "TOTAL CHARGE" in uu and "MULTIPLICITY" in uu:
                try:
                     parts = line.split()
                     self.data["charge"] = int(parts[-4])
                     self.data["mult"] = int(parts[-1])
                except: pass
            if "SCF CONVERGED" in uu or "OPTIMIZATION CONVERGED" in uu or "HURRAY" in uu:
                self.data["converged"] = True
                
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
            # Scan Step Header
            if "RELAXED SURFACE SCAN STEP" in uu_line:
                step_idx = 0
                try: step_idx = int(line.split()[-1])
                except: pass
                
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
                    elif "GEOMETRY CONVERGENCE" in uu or "CONVERGENCE CRITERIA" in uu:
                        c_idx = k + 1
                        found_any = False
                        while c_idx < next_marker and c_idx < k + 20:
                            cl = self.lines[c_idx].strip()
                            # Only break on rule if we've already found some data lines
                            if "---" in cl:
                                if found_any: break
                                else: 
                                    c_idx += 1
                                    continue
                            
                            p = cl.split()
                            if len(p) >= 4:
                                s = p[-1]
                                t = p[-2]
                                v = p[-3]
                                n = " ".join(p[:-3]).strip().lower()
                                if n and n != "item": # Skip header line
                                    conv_info[n] = {
                                        'value': v,
                                        'tolerance': t,
                                        'converged': s
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
                        'step': step_idx,
                        'energy': en,
                        'atoms': atoms,
                        'coords': coords,
                        'convergence': conv_info,
                        'gradients': step_grads
                    })
                    
            elif "OPTIMIZATION CYCLE" in uu_line:
                cycle_idx = 0
                try: 
                    parts = line.strip().split()
                    cycle_idx = int(parts[-1])
                except: pass
                
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
                    elif "GEOMETRY CONVERGENCE" in uu or "CONVERGENCE CRITERIA" in uu:
                        c_idx = k + 1
                        found_any = False
                        while c_idx < next_marker and c_idx < k + 25:
                            cl = self.lines[c_idx].strip()
                            if not cl:
                                c_idx += 1
                                continue
                            if "---" in cl:
                                if found_any: break
                                else:
                                    c_idx += 1
                                    continue
                            
                            p = cl.split()
                            if len(p) >= 4:
                                s = p[-1]
                                t = p[-2]
                                v = p[-3]
                                n = " ".join(p[:-3]).strip().lower()
                                if n and n != "item":
                                    conv_info[n] = {
                                        'value': v,
                                        'tolerance': t,
                                        'converged': s
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
                start_indices.append((i, "alpha"))
            elif "SPIN UP ORBITALS" in uu and i+1 < len(self.lines) and "---" in self.lines[i+1]:
                start_indices.append((i, "alpha"))
            elif "SPIN DOWN ORBITALS" in uu and i+1 < len(self.lines) and "---" in self.lines[i+1]:
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
                                  try:
                                      val = float(v_str)
                                      # Filter small coeffs to save space/time? 
                                      # Maybe keep all for accuracy?
                                      # Original code had > 0.01 filter.
                                      # If filter is too high? 0.01 is reasonable for isosurface.
                                      # But let's lower to 1e-4 just in case.
                                      # Or remove filter. 
                                      # Let's store all for now.
                                      self.data["mo_coeffs"][mo_idx]['coeffs'].append({
                                          "atom_idx": atom_idx,
                                          "sym": sym,
                                          "orb": orb,
                                          "coeff": val
                                      })
                                  except: pass
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
                if header.upper() in line.upper():
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
        
        # NBO Parsing
        nbo_start = -1
        # Look for the header "NATURAL ATOMIC ORBITAL AND NATURAL BOND ORBITAL ANALYSIS"
        # Or simpler "Summary of Natural Population Analysis:" or "NATURAL POPULATIONS:"
        # In ORCA 5 or NBO6 output through ORCA.
        
        # Common header in ORCA NBO output
        for i, line in enumerate(self.lines):
            if "NATURAL POPULATIONS" in line.upper() and "Summary of Natural" in self.lines[i-2 if i>2 else 0]:
                nbo_start = i
            elif "NATURAL POPULATIONS" in line.upper(): 
                 # Fallback
                 nbo_start = i

        if nbo_start != -1:
             # Parsing NBO table
             #                                     Natural Population 
             #              Natural    ---------------------------------------------
             #   Atom No    Charge        Core      Valence    Rydberg      Total 
             # ---------------------------------------------------------------------
             #     C  1    -0.12345      1.9999     4.1235     0.0100      6.1235
             # ...
             # ===============================
            
            nbo_charges = []
            curr = nbo_start + 1
            while curr < len(self.lines):
                line = self.lines[curr].strip()
                if "----------------" in line:
                    curr += 1
                    continue
                if "================" in line: break
                if not line:
                     curr += 1
                     continue
                
                parts = line.split()
                # Expected format: Atom No Charge ...
                # C 1 -0.12345 ...
                if len(parts) >= 3:
                     # Check if parst[1] is int
                     try:
                         sym = parts[0]
                         idx = int(parts[1])
                         chg = float(parts[2])
                         nbo_charges.append({
                             "atom_idx": idx - 1, # 1-based in NBO table usually
                             "atom_sym": sym,
                             "charge": chg
                         })
                     except: pass
                curr += 1
            
            if nbo_charges:
                self.data["charges"]["NBO"] = nbo_charges

        
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
        
    def parse_mo(self):
        self.data["mos"] = [] # List of dict: {id, occ, energy, spin}
        
        # Look for "ORBITAL ENERGIES"
        start_line = -1
        for i, line in enumerate(self.lines):
            if "ORBITAL ENERGIES" in line.upper():
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
        
        # Look for "ORBITAL ENERGIES" section
        start_idx = -1
        for i, line in enumerate(self.lines):
            if "ORBITAL ENERGIES" in line.upper() and i+1 < len(self.lines) and "---" in self.lines[i+1]:
                start_idx = i
                break
        
        if start_idx == -1:
            return
        
        # Skip header lines
        curr = start_idx + 3  # Skip title, separator, column headers
        
        while curr < len(self.lines):
            line = self.lines[curr].strip()
            if not line or "---" in line:
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
                    
                    self.data["orbital_energies"].append({
                        "index": orbital_idx,
                        "occupation": occupation,
                        "energy_eh": energy_eh,
                        "energy_ev": energy_ev,
                        "type": "occupied" if occupation > 0.1 else "virtual"
                    })
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
                    if "ITERATION" in uu_k and "ENERGY" in uu_k:
                        header_idx = i + k
                        break
                
                if header_idx != -1:
                    trace = []
                    idx = header_idx + 2
                    while idx < len(self.lines):
                        l_scf = self.lines[idx].strip()
                        if not l_scf or "---" in l_scf or "SUCCESS" in l_scf or "Energy Check" in l_scf:
                            break
                        
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

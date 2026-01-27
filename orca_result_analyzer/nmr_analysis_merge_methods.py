    def merge_selected_peaks(self):
        """Merge selected peaks into a single entry"""
        if len(self.selected_peak_indices) < 2:
            QMessageBox.warning(self, "Invalid Selection", "Please select at least 2 peaks to merge.")
            return
        
        # Get selected atom indices
        selected_indices = []
        for peak_idx in sorted(self.selected_peak_indices):
            if peak_idx < len(self.displayed_data):
                atom_idx = self.displayed_data[peak_idx].get('atom_idx', peak_idx)
                selected_indices.append(atom_idx)
        
        # Calculate averages
        total_sigma = 0.0
        total_delta = 0.0
        for atom_idx in selected_indices:
            # Find item in original data
            item = next((d for d in self.data if d.get('atom_idx') == atom_idx), None)
            if item:
                sigma = item.get('shielding', 0.0)
                delta = self.delta_ref + (self.sigma_ref - sigma)
                total_sigma += sigma
                total_delta += delta
        
        count = len(selected_indices)
        avg_sigma = total_sigma / count
        avg_delta = total_delta / count
        
        # Create merge group
        merge_group = {
            'indices': selected_indices,
            'avg_sigma': avg_sigma,
            'avg_delta': avg_delta
        }
        
        # Check if any indices are already merged
        for existing_group in self.merged_peaks:
            if any(idx in existing_group['indices'] for idx in selected_indices):
                reply = QMessageBox.question(
                    self, "Already Merged",
                    "Some selected peaks are already in a merged group. Replace existing merge?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply == QMessageBox.StandardButton.Yes:
                    # Remove old group
                    self.merged_peaks.remove(existing_group)
                else:
                    return
        
        # Add new merge group
        self.merged_peaks.append(merge_group)
        
        # Save to JSON
        self.save_merged_peaks()
        
        # Clear selection and refresh
        self.clear_peak_selection()
        self.recalc()
        
        QMessageBox.information(self, "Merged", f"Successfully merged {count} peaks.")
    
    def save_merged_peaks(self):
        """Save merged peaks to JSON file"""
        try:
            with open(self.merged_peaks_file, 'w') as f:
                json.dump(self.merged_peaks, f, indent=2)
        except Exception as e:
            print(f"Error saving merged peaks: {e}")
    
    def load_merged_peaks(self):
        """Load merged peaks from JSON file"""
        if os.path.exists(self.merged_peaks_file):
            try:
                with open(self.merged_peaks_file, 'r') as f:
                    self.merged_peaks = json.load(f)
            except Exception as e:
                print(f"Error loading merged peaks: {e}")
                self.merged_peaks = []
        else:
            self.merged_peaks = []

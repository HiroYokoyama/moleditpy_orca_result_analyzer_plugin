import os

PLUGIN_NAME = "ORCA Result Analyzer"
PLUGIN_VERSION = "0.3.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = "Comprehensive analyzer for ORCA output files (.out, .log). Includes Vibrational, MO, TDDFT, and NMR analysis."

from .gui import OrcaResultAnalyzerDialog
from .parser import OrcaParser

# Global reference to keep window alive
_analyzer_window = None

def initialize(context):
    """
    Initialize the ORCA Result Analyzer plugin.
    Registers file openers for .out and .log with HIGH PRIORITY (100).
    """

    def open_orca_file(path):
        global _analyzer_window
        mw = context.get_main_window()
        
        # Read file to memory
        try:
            with open(path, 'r', encoding='utf-8', errors='replace') as f:
                content = f.read()
        except Exception as e:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.critical(mw, "Error Reading File", f"Could not read file:\n{e}")
            return

        # Initialize Parser
        # Force reload to ensure latest code is used (dev mode helper)
        import importlib
        from . import parser as parser_mod
        importlib.reload(parser_mod)
        from .parser import OrcaParser
        
        parser = OrcaParser()
        parser.load_from_memory(content, path)
        
        # Close existing if open
        if _analyzer_window is not None:
            try:
                _analyzer_window.close()
            except: pass
            _analyzer_window = None
            
        # Open Dialog (Modeless)
        from .gui import OrcaResultAnalyzerDialog
        _analyzer_window = OrcaResultAnalyzerDialog(mw, parser, path, context)
        
        _analyzer_window.show()
        _analyzer_window.raise_()
        _analyzer_window.activateWindow()
        
        # Auto-load 3D structure
        _analyzer_window.load_structure_3d()

    def handle_drop(path):
        if path.lower().endswith((".out", ".log")):
            open_orca_file(path)
            return True
        return False

    # Register Opener
    # Priority 100 as requested
    context.register_file_opener(".out", open_orca_file, priority=100)
    context.register_file_opener(".log", open_orca_file, priority=100)
    
    # Register Drop Handler
    context.register_drop_handler(handle_drop, priority=100)

    # Add to Main Menu
    def menu_action():
        from PyQt6.QtWidgets import QFileDialog
        mw = context.get_main_window()
        path, _ = QFileDialog.getOpenFileName(mw, "Open ORCA Output", "", "ORCA Output (*.out *.log)")
        if path:
             open_orca_file(path)

    context.add_menu_action("Analysis/ORCA Result Analyzer", menu_action) 


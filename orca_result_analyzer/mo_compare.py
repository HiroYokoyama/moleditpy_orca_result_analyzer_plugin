"""Side-by-side comparison of up to four molecular orbitals.

Opened from the MO dialog. Each slot owns its own colours, isovalue, opacity
and style, and renders into its own namespaced pair of actors so all four
isosurfaces coexist in the host's 3D view.
"""

import json
import logging
import os

from PyQt6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)
from PyQt6.QtCore import Qt

try:
    from .vis import CubeVisualizer
except ImportError:
    try:
        from vis import CubeVisualizer
    except ImportError:
        CubeVisualizer = None

try:
    from .utils import save_json_atomic
except ImportError:
    from utils import save_json_atomic


SLOT_COUNT = 4

# Slot 1 is overwritten with the MO dialog's own colours at build time, so the
# comparison always opens showing the orbital the user was already looking at
# in the colours they chose there.
DEFAULT_COLORS = [
    ("#ff0000", "#0000ff"),
    ("#ff8c00", "#008b8b"),
    ("#9400d3", "#7fbf00"),
    ("#ff1493", "#00bfff"),
]

STYLES = ["Surface", "Wireframe", "Points"]

# Slot defaults when the MO table has no selection to inherit.
DEFAULT_TAGS = ["HOMO", "LUMO", "LUMO+1", "HOMO-1"]


def contrast_text(hex_c):
    """'black' or 'white', whichever stays readable on *hex_c*."""
    h = str(hex_c).lstrip("#")
    if len(h) != 6:
        return "black"
    try:
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    except ValueError:
        return "black"
    return "black" if (r * 299 + g * 587 + b * 114) / 1000 > 128 else "white"


class MOSlot:
    """One orbital row: its widgets and the values they hold."""

    def __init__(self, index, orbitals, colors, parent_dlg):
        self.index = index
        self.prefix = f"mo_cmp{index}"
        self.vis = None

        self.box = QGroupBox(f"Orbital {index + 1}")
        outer = QVBoxLayout(self.box)

        row1 = QHBoxLayout()
        self.check_on = QCheckBox("Show")
        row1.addWidget(self.check_on)

        self.combo_mo = QComboBox()
        for label, key, display_id, _tag in orbitals:
            self.combo_mo.addItem(label, (key, display_id))
        row1.addWidget(self.combo_mo)
        outer.addLayout(row1)

        row2 = QHBoxLayout()
        self.btn_p = QPushButton("Pos (+)")
        self.btn_n = QPushButton("Neg (-)")
        self.btn_p.clicked.connect(lambda: parent_dlg.pick_slot_color(self, "p"))
        self.btn_n.clicked.connect(lambda: parent_dlg.pick_slot_color(self, "n"))
        row2.addWidget(self.btn_p)
        row2.addWidget(self.btn_n)
        outer.addLayout(row2)

        row3 = QHBoxLayout()
        row3.addWidget(QLabel("Iso:"))
        self.spin_iso = QDoubleSpinBox()
        self.spin_iso.setRange(0.001, 1.0)
        self.spin_iso.setSingleStep(0.005)
        self.spin_iso.setDecimals(3)
        self.spin_iso.setValue(0.02)
        row3.addWidget(self.spin_iso)

        row3.addWidget(QLabel("Opacity:"))
        self.spin_opacity = QDoubleSpinBox()
        self.spin_opacity.setRange(0.0, 1.0)
        self.spin_opacity.setSingleStep(0.1)
        self.spin_opacity.setValue(0.5)
        row3.addWidget(self.spin_opacity)
        outer.addLayout(row3)

        row4 = QHBoxLayout()
        row4.addWidget(QLabel("Style:"))
        self.combo_style = QComboBox()
        self.combo_style.addItems(STYLES)
        row4.addWidget(self.combo_style)
        self.check_smooth = QCheckBox("Smooth")
        self.check_smooth.setChecked(True)
        row4.addWidget(self.check_smooth)
        outer.addLayout(row4)

        self.set_color("p", colors[0])
        self.set_color("n", colors[1])

    def set_color(self, which, hex_c):
        btn = self.btn_p if which == "p" else self.btn_n
        btn.setStyleSheet(
            f"background-color: {hex_c}; color: {contrast_text(hex_c)}; "
            "font-weight: bold;"
        )

    def color(self, which):
        btn = self.btn_p if which == "p" else self.btn_n
        style = btn.styleSheet()
        if "background-color:" in style:
            return style.split("background-color:")[1].split(";")[0].strip()
        return "#ff0000" if which == "p" else "#0000ff"

    def is_on(self):
        return self.check_on.isChecked()

    def selection(self):
        """(access key, display id) of the chosen orbital, or (None, None)."""
        data = self.combo_mo.currentData()
        if isinstance(data, (tuple, list)) and len(data) == 2:
            return data[0], data[1]
        return None, None

    def to_settings(self):
        """Appearance only. Which orbital a slot holds, and whether it is on,
        follow the loaded file's selection and are not persisted."""
        return {
            "color_p": self.color("p"),
            "color_n": self.color("n"),
            "iso": self.spin_iso.value(),
            "opacity": self.spin_opacity.value(),
            "style": self.combo_style.currentText(),
            "smooth_shading": bool(self.check_smooth.isChecked()),
        }

    def from_settings(self, data):
        """Restore appearance from *data*; unknown or malformed keys are left
        at whatever the slot was seeded with."""
        if not isinstance(data, dict):
            return
        for which, key in (("p", "color_p"), ("n", "color_n")):
            value = data.get(key)
            if isinstance(value, str) and value:
                self.set_color(which, value)
        for widget, key in ((self.spin_iso, "iso"), (self.spin_opacity, "opacity")):
            value = data.get(key)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                widget.setValue(float(value))
        style = data.get("style")
        if isinstance(style, str):
            idx = self.combo_style.findText(style)
            if idx >= 0:
                self.combo_style.setCurrentIndex(idx)
        if "smooth_shading" in data:
            self.check_smooth.setChecked(bool(data.get("smooth_shading")))


class MOCompareDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent_dlg = parent
        self.mw = getattr(parent, "mw", None)
        self.setWindowTitle("Compare Molecular Orbitals")
        self.resize(760, 520)
        self.slots = []
        self.setup_ui()

    # -- construction ------------------------------------------------------

    def collect_orbitals(self):
        """(label, access key, display id, frontier tag) per MO table row."""
        out = []
        tree = getattr(self.parent_dlg, "tree", None)
        if tree is None:
            return out
        for i in range(tree.topLevelItemCount()):
            item = tree.topLevelItem(i)
            if item is None:
                continue
            display_id = item.text(0)
            tag = item.text(1)
            label = f"{display_id}  {tag}".strip() if tag else display_id
            key = item.data(0, Qt.ItemDataRole.UserRole)
            if key is None:
                continue
            out.append((label, key, display_id, tag))
        return out

    def preselected_keys(self):
        """Orbitals highlighted in the MO table when this window opened."""
        keys = []
        tree = getattr(self.parent_dlg, "tree", None)
        if tree is None:
            return keys
        try:
            selected = tree.selectedItems() or []
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)
            return keys
        for item in selected:
            key = item.data(0, Qt.ItemDataRole.UserRole)
            if key is not None and key not in keys:
                keys.append(key)
        return keys[:SLOT_COUNT]

    def default_keys(self):
        """The frontier orbitals, in DEFAULT_TAGS order.

        Unrestricted output tags both spin channels, so a tag can repeat; the
        first row wins, which is the channel the table lists first.
        """
        by_tag = {}
        for _label, key, _display_id, tag in self.orbitals:
            if tag and tag not in by_tag:
                by_tag[tag] = key
        return [by_tag[t] for t in DEFAULT_TAGS if t in by_tag]

    def apply_selection(self, keys=None):
        """Fill the slots with *keys*, or with the table selection, or with
        the frontier orbitals — whichever is available, in that order."""
        if keys is None:
            keys = self.preselected_keys() or self.default_keys()
        positions = {key: i for i, (_l, key, _d, _t) in enumerate(self.orbitals)}

        assigned = 0
        for slot in self.slots:
            key = keys[assigned] if assigned < len(keys) else None
            pos = positions.get(key) if key is not None else None
            if pos is None:
                slot.check_on.setChecked(False)
                continue
            slot.combo_mo.setCurrentIndex(pos)
            slot.check_on.setChecked(True)
            assigned += 1

    def setup_ui(self):
        layout = QVBoxLayout(self)

        info = QLabel(
            "Pick up to four orbitals. Missing cubes are generated on Update."
        )
        layout.addWidget(info)

        self.orbitals = self.collect_orbitals()
        grid = QGridLayout()
        for i in range(SLOT_COUNT):
            colors = list(DEFAULT_COLORS[i])
            if i == 0:
                colors = [self._parent_color("p"), self._parent_color("n")]
            slot = MOSlot(i, self.orbitals, colors, self)
            if i == 0:
                self._seed_first_slot(slot)
            grid.addWidget(slot.box, i // 2, i % 2)
            self.slots.append(slot)
        layout.addLayout(grid)

        # Saved appearance wins over the seeded defaults; with nothing saved,
        # slot 1 keeps the colours it just inherited from the MO dialog.
        self.load_settings()
        # Which orbitals go where depends on the loaded file, not on saved
        # settings, so this runs last.
        self.apply_selection()

        btns = QHBoxLayout()
        self.btn_update = QPushButton("Update View")
        self.btn_update.setStyleSheet("font-weight: bold; background-color: #d0f0c0;")
        self.btn_update.clicked.connect(self.update_view)
        btns.addWidget(self.btn_update)

        btn_clear = QPushButton("Clear All")
        btn_clear.clicked.connect(self.clear_all)
        btns.addWidget(btn_clear)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btns.addWidget(btn_close)
        layout.addLayout(btns)

        self.lbl_status = QLabel("")
        layout.addWidget(self.lbl_status)

    # -- persistence -------------------------------------------------------

    def settings_path(self):
        """Shared with the MO dialog: one settings.json beside the package."""
        return os.path.join(os.path.dirname(__file__), "settings.json")

    def load_settings(self):
        path = self.settings_path()
        if not os.path.exists(path):
            return
        try:
            with open(path, "r", encoding="utf-8") as fh:
                saved = json.load(fh).get("mo_compare", {})
        except (OSError, ValueError) as e:
            logging.warning("Error loading compare settings: %s", e)
            return
        slots = saved.get("slots") if isinstance(saved, dict) else None
        if not isinstance(slots, list):
            return
        for slot, data in zip(self.slots, slots):
            slot.from_settings(data)

    def save_settings(self):
        path = self.settings_path()
        all_settings = {}
        if os.path.exists(path):
            try:
                with open(path, "r", encoding="utf-8") as fh:
                    all_settings = json.load(fh)
            except (OSError, ValueError) as _e:
                logging.warning("silenced: %s", _e)
        if not isinstance(all_settings, dict):
            all_settings = {}

        all_settings["mo_compare"] = {
            "slots": [slot.to_settings() for slot in self.slots]
        }
        try:
            save_json_atomic(path, all_settings)
        except (OSError, TypeError, ValueError) as e:
            logging.warning("Error saving compare settings: %s", e)

    def _parent_color(self, which):
        try:
            return self.parent_dlg.get_color_hex(which)
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)
            return DEFAULT_COLORS[0][0 if which == "p" else 1]

    def _seed_first_slot(self, slot):
        """Slot 1 renders with the MO dialog's own settings."""
        parent = self.parent_dlg
        try:
            slot.spin_iso.setValue(parent.spin_iso.value())
            slot.spin_opacity.setValue(parent.spin_opacity.value())
            idx = slot.combo_style.findText(parent.combo_style.currentText())
            if idx >= 0:
                slot.combo_style.setCurrentIndex(idx)
            slot.check_smooth.setChecked(parent.check_smooth.isChecked())
        except (AttributeError, RuntimeError, TypeError) as _e:
            logging.warning("silenced: %s", _e)

    # -- interaction -------------------------------------------------------

    def pick_slot_color(self, slot, which):
        from PyQt6.QtGui import QColor

        col = QColorDialog.getColor(QColor(slot.color(which)), self, "Select Color")
        if col.isValid():
            slot.set_color(which, col.name())
            self.render_all()

    def update_view(self):
        """Generate whatever cube files are missing, then draw everything."""
        missing = []
        for slot in self.slots:
            if not slot.is_on():
                continue
            key, display_id = slot.selection()
            if key is None:
                continue
            path = self._cube_path(display_id)
            if not path or not os.path.exists(path):
                if key not in missing:
                    missing.append(key)

        if not missing:
            self.render_all()
            return

        self.lbl_status.setText(f"Generating {len(missing)} cube(s)...")
        try:
            self.parent_dlg.generate_cubes(missing, on_done=self.render_all)
        except (AttributeError, RuntimeError) as e:
            logging.warning("Cube generation unavailable: %s", e)
            QMessageBox.warning(
                self, "Error", f"Could not generate the missing cubes:\n{e}"
            )
            self.lbl_status.setText("")

    def _cube_path(self, display_id):
        try:
            return self.parent_dlg.get_cube_path(display_id)
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)
            return None

    def render_all(self):
        if not CubeVisualizer:
            QMessageBox.warning(
                self,
                "Visualizer Error",
                "CubeVisualizer module not loaded.\nCheck if 'pyvista' is installed.",
            )
            return
        mw = self.mw
        if not mw:
            return

        # The MO dialog's own single-orbital actors would sit on top of the
        # comparison and double-draw slot 1.
        self._remove_actors("mo_iso")

        shown, missing = 0, 0
        for slot in self.slots:
            if not slot.is_on():
                self._remove_actors(slot.prefix)
                continue
            _key, display_id = slot.selection()
            path = self._cube_path(display_id)
            if not path or not os.path.exists(path):
                self._remove_actors(slot.prefix)
                missing += 1
                continue

            vis = CubeVisualizer(mw)
            if not vis.load_file(path):
                self._remove_actors(slot.prefix)
                continue
            vis.show_iso(
                slot.spin_iso.value(),
                opacity=slot.spin_opacity.value(),
                color_p=slot.color("p"),
                color_n=slot.color("n"),
                style=slot.combo_style.currentText().lower(),
                smooth_shading=slot.check_smooth.isChecked(),
                name_prefix=slot.prefix,
            )
            slot.vis = vis
            shown += 1

        try:
            mw.plotter.render()
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)

        msg = f"Showing {shown} orbital(s)."
        if missing:
            msg += f" {missing} cube(s) still missing."
        self.lbl_status.setText(msg)
        self.save_settings()

    def clear_all(self):
        for slot in self.slots:
            slot.check_on.setChecked(False)
            self._remove_actors(slot.prefix)
        try:
            if self.mw:
                self.mw.plotter.render()
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)
        self.lbl_status.setText("")

    def _remove_actors(self, prefix):
        try:
            plotter = self.mw.plotter if self.mw else None
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)
            return
        if plotter is None:
            return
        for suffix in ("_p", "_n"):
            try:
                plotter.remove_actor(f"{prefix}{suffix}")
            except (AttributeError, RuntimeError, KeyError) as _e:
                logging.warning("silenced: %s", _e)

    # -- teardown ----------------------------------------------------------

    def reject(self):
        # Esc must run closeEvent cleanup (QDialog.reject only hides)
        self.close()

    def closeEvent(self, event):
        self.save_settings()
        for slot in self.slots:
            self._remove_actors(slot.prefix)
        try:
            if self.mw:
                self.mw.plotter.render()
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)

        # Without this the MO dialog keeps a dead reference and reopening
        # raises the destroyed window instead of building a new one.
        try:
            self.parent_dlg.on_compare_closed()
        except (AttributeError, RuntimeError) as _e:
            logging.warning("silenced: %s", _e)

        # accept() not super().closeEvent(): QDialog.closeEvent calls reject(),
        # which is routed back through close() and would recurse.
        event.accept()

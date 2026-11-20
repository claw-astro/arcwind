"""
Figures. Every function returns a matplotlib Figure (or Axes) and accepts
``save="path.pdf"``.

=================================  ====================================================
function                           shows
=================================  ====================================================
:func:`plot_pv`                    one PV diagram (+ mask, trace, model)
:func:`plot_pv_grid`               PV diagrams of many cuts (Louvet+2018 Fig. 7)
:func:`plot_gallery`               high-contrast gallery with model and bootstrap
:func:`plot_slit_layout`           the cuts on the moment-0 / continuum map
:func:`plot_cuts`                  every fitted cut with trace and model(s)
:func:`plot_cut`                   zoom on one cut
:func:`plot_bootstrap`             bootstrap histograms
:func:`plot_profiles`              Vz, Vr, Vp, flow angle, R, Vphi vs height
:func:`plot_ring_fits`             per-cut ring fits vs height
:func:`plot_mass`                  moment map, brightness, mass and Mdot
:func:`plot_noise_check`           noise-channel diagnostics
:func:`plot_disk_rotation_check`   disk rotation sense
=================================  ====================================================
"""
from .fit import plot_bootstrap, plot_cut, plot_cuts, plot_profiles, plot_ring_fits
from .maps import plot_slit_layout
from .other import plot_disk_rotation_check, plot_mass, plot_noise_check
from .pv import plot_gallery, plot_pv, plot_pv_grid
from ._style import set_style

__all__ = ["set_style", "plot_pv", "plot_pv_grid", "plot_gallery", "plot_slit_layout", "plot_cuts",
           "plot_cut", "plot_bootstrap", "plot_profiles", "plot_ring_fits", "plot_mass",
           "plot_noise_check", "plot_disk_rotation_check"]

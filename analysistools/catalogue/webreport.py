#!/usr/bin/env python3
"""
analysistools.catalogue.webreport
------------------------------------
WebReport: a static, self-contained HTML site for sharing catalogue-derived
figures with collaborators or as paper supplementary material, without a
server and without handing out the catalogue file itself.

Uses ``plotly`` (already a core dependency, see pyproject.toml) via
``fig.to_html(full_html=False, include_plotlyjs="cdn")`` for embeddable,
interactive (pan/zoom/hover) figures. Output is a plain folder of
HTML/JS -- deployable to GitHub Pages, attached to a Zenodo deposit, or
zipped and emailed.

This is a snapshot generator, not a live dashboard: figures are rendered
once from a given ``CatalogueDataset`` (or any ``Dataset``) and the
catalogue version/provenance string used is stamped into the page footer,
so a reader always knows which release a figure came from. For an
always-current view backed by a running catalogue, see the Cowork
`create_artifact`-style approach instead -- nothing here forecloses that,
it simply isn't what a supplementary-material link needs.

See docs/master_catalogue.md ("Shareable results").
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import logging

logger = logging.getLogger(__name__)


@dataclass
class FigureSpec:
    """One figure to include in the report."""

    name: str
    title: str
    figure: Any                # a plotly.graph_objects.Figure
    description: str = ""
    group: str = "default"     # groups figures onto separate report pages


@dataclass
class WebReport:
    """Accumulates figures, then writes a small static multi-page site.

    Parameters
    ----------
    catalogue_version : str
        Stamped into every page footer for provenance.
    title : str
        Site title (e.g. "Dorcha Paper I -- Supplementary Figures").
    """

    catalogue_version: str
    title: str = "Catalogue Report"
    figures: List[FigureSpec] = field(default_factory=list)

    def add_figure(self, name: str, figure: Any, *, title: str = "",
                  description: str = "", group: str = "default"
                  ) -> "WebReport":
        self.figures.append(FigureSpec(
            name=name, title=title or name, figure=figure,
            description=description, group=group))
        return self

    def build(self, out_dir: str,
             include_plotlyjs: str = "cdn") -> Path:
        """Write ``index.html`` plus one page per figure group into
        ``out_dir``. Returns the output directory path.
        """
        raise NotImplementedError(
            "Phase 6c: for each group, embed figures via "
            "fig.to_html(full_html=False, include_plotlyjs=include_plotlyjs) "
            "into a page template with the catalogue_version/title footer; "
            "write an index.html linking every group page.")

    @classmethod
    def from_catalogue(cls, catalogue: "analysistools.api.dataset.Dataset",
                       title: str = "Catalogue Report") -> "WebReport":
        """Convenience constructor reading ``catalogue.meta`` for the
        version/provenance stamp."""
        version = catalogue.meta.get("catalogue_version", "unknown")
        return cls(catalogue_version=version, title=title)

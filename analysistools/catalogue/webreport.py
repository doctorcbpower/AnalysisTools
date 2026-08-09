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
from html import escape
from pathlib import Path
from typing import Any, Dict, List, Optional

import logging

logger = logging.getLogger(__name__)

_PAGE_TEMPLATE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body {{ font-family: -apple-system, Helvetica, Arial, sans-serif;
       margin: 0; padding: 0 2rem 3rem; color: #1a1a1a; }}
header {{ padding: 1.5rem 0 1rem; border-bottom: 1px solid #ddd; }}
h1 {{ margin: 0 0 0.25rem; }}
nav a {{ margin-right: 1rem; }}
.figure {{ margin: 2rem 0; }}
.figure h2 {{ margin-bottom: 0.25rem; }}
.figure p {{ color: #555; }}
footer {{ margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #ddd;
         color: #888; font-size: 0.85rem; }}
</style>
</head>
<body>
<header>
<h1>{title}</h1>
<nav>{nav}</nav>
</header>
{body}
<footer>{footer}</footer>
</body>
</html>
"""


def _nav_links(groups: List[str], current: Optional[str]) -> str:
    links = []
    index_label = "Index"
    index_href = "index.html"
    links.append(f'<a href="{index_href}">{index_label}</a>'
                if current is not None else f"<strong>{index_label}</strong>")
    for group in groups:
        href = f"{group}.html"
        label = escape(group)
        if group == current:
            links.append(f"<strong>{label}</strong>")
        else:
            links.append(f'<a href="{href}">{label}</a>')
    return "\n".join(links)


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

        Each group becomes its own static page (``<group>.html``) with
        every figure in that group embedded via ``fig.to_html(full_html=
        False, include_plotlyjs=...)`` -- one ``<script>`` tag per figure,
        not once per page, since ``include_plotlyjs="cdn"`` (the default)
        emits a ``<script src=...>`` reference each call; passing e.g.
        ``include_plotlyjs=False`` and hosting plotly.min.js yourself
        avoids the duplication if serving many figures matters more than
        each page being a single self-contained file. ``index.html`` links
        every group page plus a plain description list of the figures in
        it. No figures at all -- an empty report -- is still written
        (an index with no group links), not an error: a report someone
        builds and looks at, with nothing selected yet, is a normal state
        while iterating on which figures to include.
        """
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)

        by_group: Dict[str, List[FigureSpec]] = {}
        for fig_spec in self.figures:
            by_group.setdefault(fig_spec.group, []).append(fig_spec)
        groups = sorted(by_group)
        footer = f"Catalogue version {escape(self.catalogue_version)}"

        for group in groups:
            body_parts = []
            for fig_spec in by_group[group]:
                fig_html = fig_spec.figure.to_html(
                    full_html=False, include_plotlyjs=include_plotlyjs)
                desc = (f"<p>{escape(fig_spec.description)}</p>"
                       if fig_spec.description else "")
                body_parts.append(
                    f'<div class="figure">\n'
                    f"<h2>{escape(fig_spec.title)}</h2>\n{desc}\n{fig_html}\n"
                    f"</div>")
            page = _PAGE_TEMPLATE.format(
                title=escape(f"{self.title} – {group}"),
                nav=_nav_links(groups, current=group),
                body="\n".join(body_parts), footer=footer)
            (out / f"{group}.html").write_text(page)

        index_body = []
        for group in groups:
            items = "\n".join(
                f"<li>{escape(fig_spec.title)}"
                + (f" — {escape(fig_spec.description)}"
                  if fig_spec.description else "")
                + "</li>"
                for fig_spec in by_group[group])
            index_body.append(
                f'<div class="figure">\n<h2><a href="{group}.html">'
                f"{escape(group)}</a></h2>\n<ul>{items}</ul>\n</div>")
        index_page = _PAGE_TEMPLATE.format(
            title=escape(self.title), nav=_nav_links(groups, current=None),
            body="\n".join(index_body), footer=footer)
        (out / "index.html").write_text(index_page)

        logger.info("WebReport: wrote %d group page(s) + index.html to %s",
                   len(groups), out)
        return out

    @classmethod
    def from_catalogue(cls, catalogue: "analysistools.api.dataset.Dataset",
                       title: str = "Catalogue Report") -> "WebReport":
        """Convenience constructor reading ``catalogue.meta`` for the
        version/provenance stamp."""
        version = catalogue.meta.get("catalogue_version", "unknown")
        return cls(catalogue_version=version, title=title)

// QAtlas.jl documentation enhancements
//
// Three additions on every rendered page:
//
//   1. A floating "Report an issue" button fixed at the top-right of
//      the viewport, always visible, linking to a prefilled GitHub issue.
//
//   2. A banner at the top of the article stating that the documentation
//      is an AI-assisted draft.
//
//   3. A small "report" link appended to every h2 / h3 heading that
//      opens a GitHub issue prefilled with the section title and
//      page URL.
//
// All hooks target the GitHub repository `sotashimozono/QAtlas.jl`.

(function () {
    const REPO = "sotashimozono/QAtlas.jl";

    function issueUrl(title, body) {
        const params = new URLSearchParams({
            title: title,
            body: body,
            labels: "docs",
        });
        return "https://github.com/" + REPO + "/issues/new?" + params.toString();
    }

    function init() {
        const article = document.querySelector("article.docs-content") ||
                        document.querySelector("article");
        if (!article) return;

        const pageUrl = window.location.href;

        // ── 1. Fixed top-right "Report an issue" button ─────────────────
        const fab = document.createElement("a");
        fab.href = issueUrl(
            "[docs] error report",
            "Page: " + pageUrl + "\n\nIssue:\n"
        );
        fab.target = "_blank";
        fab.rel = "noopener";
        fab.className = "report-fab";
        fab.title = "Report an issue with this page";
        fab.setAttribute("aria-label", "Report an issue with this page");
        fab.textContent = "Report an issue";
        document.body.appendChild(fab);

        // ── 2. Top-of-page banner ────────────────────────────────────────
        const banner = document.createElement("div");
        banner.className = "ai-draft-banner";
        banner.innerHTML =
            '<strong>AI-assisted draft.</strong> ' +
            'This documentation is largely generated with LLM assistance ' +
            '(Claude) and every derivation is still being independently ' +
            'reviewed. If you spot an error, please ' +
            '<a href="' + issueUrl(
                "[docs] error report",
                "Page: " + pageUrl + "\n\nIssue:\n"
            ) + '" target="_blank" rel="noopener">open an issue</a> ' +
            'or submit a pull request via the ' +
            '<em>Edit on GitHub</em> link at the bottom of the page. ' +
            'Feedback and corrections are very welcome.';
        article.insertBefore(banner, article.firstChild);

        // ── 3. Per-section report links ──────────────────────────────────
        article.querySelectorAll("h2, h3").forEach(function (h) {
            if (h.closest(".ai-draft-banner")) return;

            const sectionTitle = h.textContent.trim().replace(/\s+/g, " ");
            if (!sectionTitle) return;

            const link = document.createElement("a");
            link.href = issueUrl(
                "[docs] " + sectionTitle,
                "Section: **" + sectionTitle + "**\n" +
                "Page: " + pageUrl + "\n\nIssue:\n"
            );
            link.target = "_blank";
            link.rel = "noopener";
            link.className = "report-section-link";
            link.title = "Report an issue with this section";
            link.setAttribute("aria-label", "Report an issue with this section");
            link.textContent = "report";
            h.appendChild(link);
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();

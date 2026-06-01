// QAtlas.jl documentation enhancements
//
// Two additions on every rendered page:
//
//   1. A banner at the top stating that the documentation is an
//      AI-assisted draft and that error reports / PRs are welcome,
//      with a direct link to open a new issue prefilled with the
//      page URL.
//
//   2. A small "report" link appended to every h2 / h3 heading that
//      opens a GitHub issue prefilled with the section title and
//      page URL, so readers who spot something wrong in a single
//      section can report it with one click.
//
// Both hooks target the GitHub repository `sotashimozono/QAtlas.jl`.

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

        // ── 1. Top-of-page banner ────────────────────────────────────────
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

        // ── 2. Per-section report links ──────────────────────────────────
        article.querySelectorAll("h2, h3").forEach(function (h) {
            // Skip the banner itself and Documenter-generated headings
            // that shouldn't have a report link.
            if (h.closest(".ai-draft-banner")) return;

            // Heading text may contain math, anchors, etc.  Use the
            // textContent as the issue title, but keep it short.
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

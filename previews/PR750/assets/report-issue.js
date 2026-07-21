(function () {
    const REPO = "QAtlasHub/QAtlas.jl";

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

        // ── AI-draft banner ──────────────────────────────────────────────
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
            '<em>Edit on GitHub</em> link at the bottom of the page.';
        article.insertBefore(banner, article.firstChild);

        // ── Per-section right-aligned "report" text ──────────────────────
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
            link.setAttribute("aria-label", "Report an issue with this section");
            link.textContent = "report";
            h.insertBefore(link, h.firstChild);
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();

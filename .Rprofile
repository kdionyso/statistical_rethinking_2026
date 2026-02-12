if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
    if ("httpgd" %in% .packages(all.available = TRUE)) {
        options(vsc.plot = FALSE)
        options(device = function(...) {
            httpgd::hgd(silent = TRUE)
            .vsc.browser(httpgd::hgd_url(), viewer = "Beside")
        })
    }
}

# if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
#     if (requireNamespace("httpgd", quietly = TRUE)) {
#         options(vsc.plot = FALSE)
#         options(device = function(...) {
#             httpgd::hgd(silent = TRUE)
#             .vsc.browser(httpgd::hgd_url(history = FALSE), viewer = "Beside")
#         })
#     }
# }

# if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
#     if (requireNamespace("httpgd", quietly = TRUE)) {
#         # options(vsc.plot = FALSE)
#         options(vsc.plot = FALSE) # , vsc.use_httpgd = FALSE)
#         options(device = function(...) {
#             httpgd::hgd(silent = TRUE)
#             .vsc.browser(httpgd::hgd_url(history = TRUE), viewer = "Beside")
#         })
#         # options(device = function(...) {
#         #     httpgd::hgd(silent = TRUE)
#         #     .vsc.browser(httpgd::hgd_url(history = FALSE), viewer = "Beside")
#         # })
#     }
# }

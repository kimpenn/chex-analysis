## Author: Youtao Lu <luyoutao@sas.upenn.edu>
 
## Copyright (c) 2017-2023, Youtao Lu and Junhyong Kim, Department of Biology, University of Pennsylvania
## Copyright (c) 2017-2023, James Eberwine, Perelman School of Medicine, University of Pennsylvania
## All Rights Reserved.

## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
#######################################################################
if (!exists("Tools") || is.environment(Tools)) { 
    Tools <- new.env(parent = emptyenv())
}

local({
    .VERSION = "0.3"
    library("tibble")
    options(tibble.width = Inf)
    library("scales") 
    library("stringr")
    library("openxlsx")

    is_online <- function(host = "\\.upenn\\.edu$") {
        grepl(host, system("hostname", intern = TRUE), ignore.case = TRUE)
    }


    list2data.frame <- function(list, fill = "") {
        n <- max(sapply(list, length))
        as.data.frame(sapply(list, function(l) c(l, rep(fill, n - length(l))), simplify = TRUE), 
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    write_xlsx <- function(dflist, filename, sheetnames = NULL, row.names = TRUE) {
        wb <- openxlsx::createWorkbook()
        k <- 0
        if (is.null(sheetnames)) {
            if (is.null(names(dflist))) { stop("No names in data.frame list nor sheetnames provided!") }
            sheetnames <- names(dflist)
        }
        for (x in sheetnames) {
            k <- k + 1
            openxlsx::addWorksheet(wb, sheetName = substr(x, 1, 31))
            df <- dflist[[k]]
            if (!is.data.frame(df)) { df <- as.data.frame(df) }
            openxlsx::writeDataTable(wb, sheet = k, x = df, keepNA = TRUE, headerStyle = NULL, tableStyle = "none", bandedRows = FALSE, bandedCols = FALSE, withFilter = TRUE, rowNames = row.names)
        }
        openxlsx::saveWorkbook(wb, file = filename, overwrite = TRUE)
    }

    read_xlsx <- function(filename, sheetnames, ...) {
        dfs <- lapply(sheetnames, function(sheet) {
        openxlsx::read.xlsx(xlsxFile = filename, sheet = sheet, ...)
        })
        names(dfs) <- sheetnames
        dfs
    }

    for (.obj in ls()) {
        assign(.obj, get(.obj), envir = Tools)
    }
})

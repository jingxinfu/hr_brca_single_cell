options(warn = 2)

installed_p <- installed.packages()[, "Package"]

#-------------------- Packages from R CRAN --------------------
r_packages <- c(
                "IRkernel"
)
for (p in r_packages) {
    if (!(p %in% installed_p)) {
        tryCatch(
            expr = {
                install.packages(p, dependencies = TRUE, repos = "http://cran.rstudio.com/", quiet = FALSE)
                if(p=='IRkernel'){
                    IRkernel::installspec(displayname = 'scpipe_v4 R 4.1',user=TRUE)
                }
            },
            error = function(e) {
                print(e)
                quit(status = 1)
            }
        )
    }
}
#-------------------- Packages from Remotes --------------------
remotes_packages <- c(
)

for (p in remotes_packages) {
    p_name <- strsplit(p, "\\/")[[1]][1]
    version <- strsplit(p, "\\/")[[1]][2]
    if (!(p %in% installed_p)) {
        tryCatch(
            expr = {
                remotes::install_version(p_name, version = version, repos = "http://cran.rstudio.com/", quiet = FALSE)
            },
            error = function(e) {
                print(e)
                quit(status = 1)
            }
        )
    }
}

#-------------------- Packages from Github --------------------
github_packages <- c()

for (p in github_packages) {
    p_name <- strsplit(p, "\\/")[[1]][2]
    print(p)
    if (!(p_name %in% installed_p)) {
        tryCatch(
            expr = {
                devtools::install_github(p, quiet = FALSE)
            },
            error = function(e) {
                print(e)
                quit(status = 1)
            }
        )
    }
}

#-------------------- Packages from Bioconductor --------------------
bio_packages <- c(
)

for (p in bio_packages) {
    if (!(p %in% installed_p)) {
        tryCatch(
            expr = {
                BiocManager::install(p, quiet = FALSE)
            },
            error = function(e) {
                print(e)
                quit(status = 1)
            }
        )
    }
}

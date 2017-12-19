dist.dna_schultz <- function (x, model = "K80", variance = FALSE, gamma = FALSE, 
    pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE) 
{
    MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", 
        "TN93", "GG95", "LOGDET", "BH87", "PARALIN", "N", "TS", 
        "TV", "INDEL", "INDELBLOCK")
    #N is model 13
    imod <- pmatch(toupper(model), MODELS)
    if (is.na(imod)) 
        stop(paste("'model' must be one of:", paste("\"", MODELS, 
            "\"", sep = "", collapse = " ")))
    if (imod == 11 && variance) {
        warning("computing variance not available for model BH87")
        variance <- FALSE
    }
    if (gamma && imod %in% c(1, 5:7, 9:17)) {
        warning(paste("gamma-correction not available for model", 
            model))
        gamma <- FALSE
    }
    if (is.list(x)) 
        x <- as.matrix(x)
    nms <- dimnames(x)[[1]]
    n <- dim(x)
    s <- n[2]
    n <- n[1]
        if (s * n > 2^31 - 1) 
        stop("dist.dna() cannot handle more than 2^31 - 1 bases")
    if (imod %in% c(4, 6:8)) {
        BF <- if (is.null(base.freq)) 
            base.freq(x)
        else base.freq
     }
    else BF <- 0
    if (imod %in% 16:17) 
        pairwise.deletion <- TRUE
    if (!pairwise.deletion) {
        keep <- .C(GlobalDeletionDNA, x, n, s, rep(1L, s))[[4]]
        x <- x[, as.logical(keep)]
        s <- dim(x)[2]
    }

    Ndist <- if (imod == 11) 
        n * n
    else n * (n - 1)/2
    var <- if (variance) 
        double(Ndist)
    else 0
    if (!gamma) 
        gamma <- alpha <- 0
    else {
        alpha <- gamma
        gamma <- 1
    }
    
    print(paste(as.integer(n), as.integer(s)))
    d <- .C(dist_dna, x, as.integer(n), as.integer(s), imod, 
        double(Ndist), BF, as.integer(pairwise.deletion), as.integer(variance), 
        var, as.integer(gamma), as.double(alpha), NAOK = TRUE)
    if (variance) 
        var <- d[[9]]
    d <- d[[5]]
    if (imod == 11) {
        dim(d) <- c(n, n)
        dimnames(d) <- list(nms, nms)
    }
    else {
        attr(d, "Size") <- n
        attr(d, "Labels") <- nms
        attr(d, "Diag") <- attr(d, "Upper") <- FALSE
        attr(d, "call") <- match.call()
        attr(d, "method") <- model
        class(d) <- "dist"
        if (as.matrix) 
            d <- as.matrix(d)
    }
    if (variance) 
        attr(d, "variance") <- var
    d
}
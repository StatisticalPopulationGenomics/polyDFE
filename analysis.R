source("../../../../../../GitHub/polyDFE/postprocessing.R")

###############################################################################
# investigate polyDFE output
###############################################################################
# parse polyDFE output
est = parseOutput("central_chimp_A.txt")
# est has 4 different runs of polyDFE
length(est)
# each run contains information about 
# the input file that polyDFE was ran on
# the DFE model used, the best likelihood found and the corresponding gradient
# the values of all parameters (including those that were fixed and not estimated)
# the parameters that where estimated, the number of sampled sequences in the data n
# various expectations under the parameters and the estimated alpha
names(est[[1]])
# the input file used
est[[1]]$input
# the model used
est[[1]]$model
# the likelihood and gradient
est[[1]][c("lk", "grad")]
# the values of all parameters
est[[1]]$values
# the parameters polyDFE estimated
est[[1]]$estimated
# the values of the estimated parameters
est[[1]]$values[est[[1]]$estimated]
# est contains, for both neutral and selected sites,
rownames(est[[1]]$expec)
# the expected SFS, divergence counts (if relevant) and misattributed polymorphism
colnames(est[[1]]$expec)
# the estimated alpha
est[[1]]$alpha

# what are the estimated gradients
est = c(parseOutput("central_chimp_A.txt"), parseOutput("central_chimp_C.txt"),
        parseOutput("central_chimp_Del.txt"))
# what are the gradients?
grad = sapply(est, function(e) e$grad)
names(grad) = sapply(est, getModelName)
grad
# which gradients are larger than 0.001?
grad[which(grad > 0.001)]


###############################################################################
# discretize the DFE 
###############################################################################
getDiscretizedDFE(est[[1]])
getDiscretizedDFE(est[[1]], c(-100, -10, -1, 0, 1))

# we can easily calculate the discretized DFE for different models
dfe = t(sapply(est, getDiscretizedDFE, c(-100, -10, -1, 0, 1)))
# give names to the rows in dfe to use in the barplot legend
rownames(dfe) = sapply(est, getModelName)
dfe
# initialize plotting area
pdf("fig1.pdf", width = 12.6, height = 5)
par(mfrow = c(1, 2), mar = c(2, 2, 4, 0) + 0.1, xpd = TRUE)
barplot(dfe, beside = TRUE, legend.text = TRUE, 
        angle = rep(c(0, 45, -45), each = 4), density = 20, col = 1:4)
# applying log scale on the y-axis enables an easier visual comparison
# for that, we need to replace 0s with NAs
dfe[which(dfe == 0)] = NA
col = rainbow(n = nrow(dfe), end = 0.9)
barplot(dfe, beside = TRUE, log = 'y', col = col)
# draw x axis
axis(1, lwd.tick = 0, labels = FALSE)
# add legend
legend(5, 9, legend = rownames(dfe), fill = col, bty = "n", ncol = 3)
dev.off()

###############################################################################
# calculate alpha
###############################################################################
est = parseOutput("central_chimp_C.txt")
# extract alpha for all the models in est
sapply(est, function(e) e$alpha)

# when a deleterious DFE is inferred, alpha_dfe is always 0
est = parseOutput("central_chimp_Del.txt")
sapply(est, function(e) e$alpha["alpha_dfe"])

# calculate both alpha_div and alpha_dfe using R
est = c(parseOutput("central_chimp_A.txt"), parseOutput("central_chimp_C.txt"),
        parseOutput("central_chimp_Del.txt"))
# first, parse divergence data from the polyDFE input file
div = parseDivergenceData("central_chimp_sfs")
# extract and calculate alpha for all the models in est
alpha = sapply(est, function(e) c(e$alpha["alpha_dfe"], 
                                  "R alpha_dfe" = estimateAlpha(e), 
                                  e$alpha["alpha_div"], 
                                  "R alpha_div" = estimateAlpha(e, div = div)))
summary(alpha["alpha_dfe", ] - alpha["R alpha_dfe", ])
summary(alpha["alpha_div", ] - alpha["R alpha_div", ])

# control the correction for misattributed polymorphism
# calculate alpha with and without the correction for all the models in est
alpha = sapply(est, function(e) c("+ corr" = estimateAlpha(e, div = div), 
                                  "- corr" = estimateAlpha(e, div = div, poly = FALSE)))
summary(alpha["+ corr", ] - alpha["- corr", ])

# control the S_sup limits
est = c(parseOutput("central_chimp_A.txt"), parseOutput("central_chimp_C.txt"))
alpha = sapply(est, function(e) c("supLimit = 0" = estimateAlpha(e), 
                                  "supLimit = 5" = estimateAlpha(e, supLimit = 5), 
                                  "supLimit = 10" = estimateAlpha(e, supLimit = 10)))
summary(alpha["supLimit = 0", ] - alpha["supLimit = 5", ])
summary(alpha["supLimit = 5", ] - alpha["supLimit = 10", ])

###############################################################################
# hypothesis testing
###############################################################################
# compare sequentially models in two different files
compareModels("central_chimp_A.txt", "central_chimp_Del.txt")

# compare models within the same file
est = parseOutput("central_chimp_A.txt")
compareModels(est[3], est[4])$LRT

# enforce nestedness
compareModels("central_chimp_C.txt", "central_chimp_Del.txt")$LRT
compareModels("central_chimp_C.txt", "central_chimp_Del.txt", nested = TRUE)

###############################################################################
# r can account for polarization error
###############################################################################
# when r is not estimated, eps_an is estimated to be non-zero
# when r is estimated, eps_an is estimated to be zero
# AIC and LRT support this
aic = list()
for (f in c("central_chimp_A.txt", "central_chimp_C.txt", "central_chimp_Del.txt")) 
{
  est = parseOutput(f)
  without_r = compareModels(est[1], est[2])
  with_r = compareModels(est[3], est[4])
  aic[[f]] = rbind("- r" = c(without_r$AIC[, "AIC model 1"], without_r$AIC[, "AIC model 2"], 
                             without_r$LRT[, "p-value"], est[[2]]$values["eps_an"]),
                   "+ r" = c(with_r$AIC[, "AIC model 1"], with_r$AIC[, "AIC model 2"], 
                             with_r$LRT[, "p-value"], est[[4]]$values["eps_an"]))
}
aic

###############################################################################
# model averaging
###############################################################################
est = c(parseOutput("central_chimp_A.txt"), 
        parseOutput("central_chimp_C.txt"), 
        parseOutput("central_chimp_Del.txt"))
# obtain the AIC weights for all the models in est
aic = getAICweights(est)
aic
# calculate model-averaged alpha
alpha_div = sum(sapply(1:length(est), 
                       function(i) aic[i, "weight"] * est[[i]]$alpha["alpha_div"]))
alpha_dfe = sum(sapply(1:length(est), 
                       function(i) aic[i, "weight"] * est[[i]]$alpha["alpha_dfe"]))
alpha_div
alpha_dfe

###############################################################################
# should divergence data be used?
###############################################################################
# delta AIC and alpha values are impacted by type of data used
est = list(c(parseOutput("central_chimp_A.txt"), 
             parseOutput("central_chimp_C.txt"),
             parseOutput("central_chimp_Del.txt")),
           c(parseOutput("central_chimp_A_no_div.txt"),
             parseOutput("central_chimp_C_no_div.txt"),
             parseOutput("central_chimp_Del_no_div.txt")))
for (i in 1:2)
{
    # calculate aic weights
    aic = getAICweights(est[[i]])
    # set model names to the rows
    rownames(aic) = sapply(est[[i]], getModelName)
    aic = cbind(aic, t(sapply(est[[i]], function(e) e$alpha)))
    # re-order after delta aic
    aic = aic[order(aic[, "delta AIC"]), ]
    # visualize aic and estimated alpha
    print(aic)
    # calculate model-averaged alpha_dfe
    print(sum(apply(aic, 1, function(a) a["weight"] * a["alpha_dfe"])))
}

# compare visually the SFS
# parse the observed SFS
sfs = parseSFSData("central_chimp_sfs")
sfs
# parse expected SFS from polyDFE for best models according to AIC
est = parseOutput("central_chimp_C.txt")[[3]]   
sfs_div = est$expec
sfs_no_div = parseOutput("central_chimp_A_no_div.txt")[[3]]$expec
# normalize expected SFS by multiplying with the number of sites analyzed
sfs_div = sfs_div * sfs[, "length_sfs"]
sfs_no_div = sfs_no_div * sfs[, "length_sfs"]
# sample size
n = max(grep("E[P", colnames(sfs_div), fixed = TRUE))
# calculate range of the y-axis for the plot
m = min(sfs[, 1:n], sfs_div[, 1:n], sfs_no_div[, 1:n])
M = max(sfs[, 1:n], sfs_div[, 1:n], sfs_no_div[, 1:n])
# initialize plot, make it log-scale on the y-axis
pdf("fig2.pdf", width = 12.6, height = 4)
par(mfrow = c(1, 1), xpd = FALSE, mar = c(4, 4, 0, 0) + 0.1)
plot(1, 1, type = "n", ylim = c(m, M), xlim = c(1, n), 
     xlab = "i", ylab = "SFS", log = "y")
# create background colors to easily differentiate between entries in the SFS
xleft = seq(from = 1, to = n, by = 2)
rect(xleft = xleft - 0.5, ybottom = m * 0.001, 
     xright = xleft + 0.5, ytop = M * 100,
     col = gray(0.9), border = NA)
# re-draw box around plot
box()
# plot both the neutral and selected SFS
for (i in 1:2)
{
  # first, the expected SFS when divergence data was used
  points(1:n - 0.15, sfs_div[i, 1:n], col = "red", pch = 19 + i)
  # then the expected SFS when divergence data was not used
  points(1:n + 0.15, sfs_no_div[i, 1:n], col = "blue", pch = 19 + i)
  # then the observed SFS
  segments(x0 = 1:n - 0.35, x1 = 1:n + 0.35, y0 = sfs[i, 1:n], col = "black")
}
# add legend
legend("topright", bty = "n", pch = c(NA, 21, 21, 20, 21), lwd = c(1, NA, NA, NA, NA),
       col = c("black", "red", "blue", "black", "black"), 
       legend = c("obs", "expec with div", "expec without div", "neutral", "selected"))
dev.off()

# perfrom goodness-of-fit test
sfs_both = rbind(sfs_div, sfs_no_div)
fit = sapply(1:nrow(sfs_both), 
             function(i) 
                 unlist(chisq.test(x = sfs_both[i, 1:n], 
                                   p = sfs[(i + 1) %% 2 + 1, 1:n], 
                                   rescale.p = TRUE)[c("statistic", "p.value")]))
colnames(fit) = paste(rep(c("div", "no div"), each = 2), rownames(sfs_div))
fit["p.value", ]
# check overall fit
c(div = sum(fit["statistic.X-squared", 1:2]), 
  no_div = sum(fit["statistic.X-squared", 3:4]))


# use model-averaged expected SFS
est = list(est_div, est_no_div)
expec_sfs = list()
for (i in 1:2)
{
    # calculate aic weights
    aic = getAICweights(est[[i]])
    # calculate model-averaged sfs
    expec = lapply(1:length(est[[i]]), 
                   function(j) aic[j, "weight"] * est[[i]][[j]]$expec)
    expec_sfs[[i]] = sapply(1:n, 
                            function(i) 
                                rowSums(sapply(expec, function(l) l[, i])))
    # normalize expected SFS by multiplying with the number of sites analyzed
    expec_sfs[[i]] = expec_sfs[[i]] * sfs[, "length_sfs"]
}

# perfrom goodness-of-fit test
sfs_both = rbind(sfs_div, sfs_no_div)
fit = sapply(1:nrow(sfs_both), 
             function(i) 
                 unlist(chisq.test(x = sfs_both[i, 1:n], 
                                   p = sfs[(i + 1) %% 2 + 1, 1:n], 
                                   rescale.p = TRUE)[c("statistic", "p.value")]))
colnames(fit) = paste(rep(c("div", "no div"), each = 2), rownames(sfs_div))
fit["p.value", ]
# check overall fit
c(div = sum(fit["statistic.X-squared", 1:2]), 
  no_div = sum(fit["statistic.X-squared", 3:4]))

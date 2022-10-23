###### flux analysis: functions involving running simple optimization problems ######


get.opt.flux <- function(model, rxns="biomass", coefs=1, dir="max", ko=NULL, xopt=FALSE, solv.pars=get.pars("lp", list())) {
  # get the max or min flux of a single reaction in the model, given in rxns as indices or IDs as in model$rxns; if rxns="biomass" or is a single string containing the substing "biomass" (the default), optimize the biomass reaction as found by get.biomass.idx using rxns as the regex
  # or, get the max or min flux of a linear combination of reactions, with the corresponding linear coefficients given in coefs
  # to knockout reaction(s), pass KO as their indices or IDs; the lb's and ub's of these reactions will be set to 0
  # if xopt=TRUE, the optimal xopt vector will also be returned, then will return a list of list(obj, xopt); otherwise just the optimal obj value
  # solv.pars: a list of control parameters passed to solver, the default is to use the default LP pars

  c <- rep(0, ncol(model$S))
  if (length(rxns)==1 && grepl("biomass", rxns, ignore.case=TRUE)) {
    i <- get.biomass.idx(model, rgx=rxns)
  } else i <- all2idx(model, rxns)
  c[i] <- coefs
  lb <- model$lb
  ub <- model$ub
  if (!is.null(ko)) {
    ko <- all2idx(model, ko)
    lb[ko] <- 0
    ub[ko] <- 0
  }

  res <- solve.model(model, csense=dir, c=c, lb=lb, ub=ub, pars=solv.pars)[[1]]

  if (xopt) {
    return(list(obj=res$obj, xopt=res$xopt))
  } else return(res$obj)
}

fba <- get.opt.flux # flux balance analysis as a synonym of get.opt.flux

get.opt.fluxes <- function(model, rxns="all", coefs=1, dir="max", nc=1L, solv.pars=get.pars("lp", list())) {
  # a wrapper around get.opt.flux for getting the min or max fluxes of multiple reactions (or multiple linear combinations of reactions) at once (if rxns and coefs provide as lists, each element corresponds to a linear combination)
  # return a named vector with names being the corresponding rxn ID as in model$rxns, or names(rxns) if rxns is provided as a list

  if (!is.list(rxns)) {
    if (length(rxns)==1 && rxns=="all") {
      rxns <- 1:length(model$rxns) # I use model$rxns instead of ncol(S) since S can contain extra columns
    } else rxns <- all2idx(model, rxns)
    names(rxns) <- model$rxns[rxns]
    if (length(coefs)==1) coefs <- rep(coefs, length(rxns))
  } else if (is.list(coefs)) {
    rxns <- lapply(rxns, all2idx, model=model)
  } else {
    stop("If `rxns` is provided as a list, `coefs` should be a matched list as well.")
  }

  res <- unlist(pbmcapply::pbmclapply(1:length(rxns), function(i) get.opt.flux(model=model, rxns=rxns[[i]], coefs=coefs[[i]], dir=dir, solv.pars=solv.pars), mc.cores=nc))
  names(res) <- names(rxns)
  res
}

fva <- function(model, rxns="all", coefs=1, nc=1L, biomass=NULL, biomass.rgx="biomass", solv.pars=get.pars("lp", list())) {
  # flux variability analysis, for rxns given as indices of IDs as in model$rxns; "all" for all rxns
  # alternatively, can specify linear combinations of rxns by passing lists to rxns and coefs
  # biomass: whether to require minimal biomass; NULL for no constraint, or a number in [0,1] meaning requiring at least this fraction of max biomass
  # biomass.rgx: regex used to find the biomass reaction
  # return a named matrix with two columns "min" and "max", rxns in the rows

  if (!is.list(rxns)) {
    if (length(rxns)==1 && rxns=="all") {
      rxns <- 1:length(model$rxns) # I use model$rxns instead of ncol(S) since S can contain extra columns
    } else rxns <- all2idx(model, rxns)
  }

  if (!is.null(biomass)) {
    if (!is.numeric(biomass) || biomass<0 || biomass>1) stop("Invalid biomass requirement.")
    bm.idx <- get.biomass.idx(model, biomass.rgx)
    model <- set.rxn.bounds(model, bm.idx, lbs=biomass, relative=TRUE, nc=1L, solv.pars=solv.pars)
  }
  message("FVA: computing minimal fluxes, progress:")
  min <- get.opt.fluxes(model, rxns, coefs, "min", nc, solv.pars)
  message("FVA: computing maximal fluxes, progress:")
  max <- get.opt.fluxes(model, rxns, coefs, "max", nc, solv.pars)
  message("Done FVA.")

  if (is.list(rxns)) res <- data.table(id=names(rxns), vmin=min, vmax=max) else res <- data.table(id=rxns, rxn=model$rxns[rxns], vmin=min, vmax=max)
  res
}

fva1 <- function(model, rxns="all", coefs=1, nc=1L, gap=NULL, agap=NULL, keep.solv.out=FALSE, solv.pars=get.pars("lp", list())) {
  # flux variability analysis after solving model, for rxns given as indices of IDs as in model$rxns; "all" for all rxns; or specify combinations of rxns by passing lists to rxns and coefs
  # model is an arbitrary model, will first solve it to get the optimal objective value, then do FVA while requiring the objective function is at the optimal value, subject to a small margin as specified in gap and agap
  # gap and agap: relative and absolute margin for the optimal objective value constraint, i.e. |obj-obj.opt|<=agap AND |obj-obj.opt|<=|gap*obj.opt|, if one is NULL then that one is not used, default to both NULL means gap==agap==0
  # return same as fva(); but if keep.solv.out=TRUE, return list(fva.res, solver.out)

  # solve the model to get optimal objective value
  if (.pkg.var$solver=="rcplex") solv.pars$nsol <- 1 else if (.pkg.var$solver=="gurobi") solv.pars$PoolSolutions <- 1
  solv.out <- solve.model(model, pars=solv.pars)
  if (solv.out[[1]]$stat %in% .pkg.const$infeas.stat) stop("Stopped due to infeasible solution.")
  obj.opt <- solv.out[[1]]$obj
  if (!is.null(gap) && is.null(agap)) {
    lb <- obj.opt - abs(obj.opt*gap)
    ub <- obj.opt + abs(obj.opt*gap)
  } else if (is.null(gap) && !is.null(agap)) {
    lb <- obj.opt - agap
    ub <- obj.opt + agap
  } else if (!is.null(gap) && !is.null(agap)) {
    lb <- obj.opt - min(agap, abs(obj.opt*gap))
    ub <- obj.opt + min(agap, abs(obj.opt*gap))
  } else lb <- ub <- obj.opt
  model.opt <- add.constraint(model, 1:length(model$c), model$c, lb, ub)

  # forcing optimal objective value, do fva
  if (.pkg.var$solver=="rcplex") solv.pars$trace <- 0 else if (.pkg.var$solver=="gurobi") solv.pars$OutputFlag <- 0
  res <- fva(model.opt, rxns=rxns, coefs=coefs, nc=nc, solv.pars=solv.pars)

  if (keep.solv.out) res <- list(fva.res=res, solver.out=solv.out[[1]])
  res
}

get.essential.rxns <- function(model, bm.lb.rel=0.1, nc=1L, bm.rgx="biomass", solv.pars=get.pars("lp", list())) {
  # get the indices of essential reactions: those whose KO reduces biomass by > 1-bm.lb.rel (i.e. biomass_KO < biomass_WT * bm.lb.rel)

  bm0 <- get.opt.flux(model, bm.rgx, solv.pars=solv.pars)
  bms <- unlist(parallel::mclapply(1:length(model$rxns), function(i) get.opt.flux(model, bm.rgx, ko=i, solv.pars=solv.pars), mc.cores=nc))
  ess.idx <- which(bms<bm.lb.rel*bm0)
}

run.ko.screen <- function(model, rxns="all+ctrl", f, ..., nc=1L, simplify=TRUE) {
  # a wrapper function to run algorithm `f` under the original wildtype model as well as models where each rxn in rxns is knocked out (i.e. both lb and ub set to 0)
  # model is the formulated model for running the algorithm `f`; it should be derived from the original metabolic model by adding columns and rows to the S matrix, such that indices of the rxns and mets do not change
  # rxns are the indices or IDs or rxns to include in the KO screen, by default all rxns plus the control wild-type model; rxns="all" to run for all rxns w/o the ctrl; if using indices, 0 means ctrl; if using IDs, "ctrl", means ctrl; or a list, each element containing multiple reactions to KO at the same time
  # ... are further parameters passed to f()
  # nc: the number of cores for running across rxns
  # simplify: if FALSE, return a list from mclapply across rxns; if TRUE, then assume the run for each rxn yields a summary data.table, and will try to combine them with rbindlist

  if (is.list(rxns)) {
    if (is.null(names(rxns))) names(rxns) <- 1:length(rxns)
    rxns <- lapply(rxns, function(x) {
      if (length(x==1) && x=="ctrl") 0 else all2idx(model, x)
    })
  } else if (length(rxns)==1 && rxns=="all+ctrl") {
    rxns <- 0:length(model$rxns) # I use model$rxns instead of ncol(S) since S can contain extra columns
    names(rxns) <- c("ctrl", model$rxns)
  } else if (length(rxns)==1 && rxns=="all") {
    rxns <- 1:length(model$rxns)
    names(rxns) <- model$rxns
  } else if (length(rxns)==1 && rxns=="ctrl") {
    rxns <- 0
    names(rxns) <- "ctrl"
  } else if (is.character(rxns)) {
    tmp <- rxns
    rxns <- match(rxns, model$rxns)
    names(rxns) <- tmp
    if ("ctrl" %in% names(rxns)) rxns["ctrl"] <- 0
  } else if (is.numeric(rxns)) {
    names(rxns)[rxns!=0] <- model$rxns[rxns[rxns!=0]]
    names(rxns)[rxns==0] <- "ctrl"
  }

  tmp <- 1:length(rxns)
  names(tmp) <- names(rxns)
  message("Begin KO screen, progress:")
  res <- pbmcapply::pbmclapply(tmp, function(i) {
    m <- model
    m$lb[rxns[[i]]] <- 0 # if rxns[[i]]==0, lb will not be changed (i.e. the control)
    m$ub[rxns[[i]]] <- 0 # if rxns[[i]]==0, ub will not be changed (i.e. the control)
    f(m, ...)
  }, mc.cores=nc)

  if (simplify) {
    if (any(!sapply(res, is.data.table))) {
      warning("Outputs contain non-data.table objects, cannot simplify.", immediate.=TRUE)
    } else res <- cbind(data.table(id=rxns), rbindlist(res, idcol="rxn"))
  }
  res
}


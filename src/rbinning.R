# binning functions
# Deloitte 2018

preprocess = function(x) {
  classx = class(x)[1]
  if (classx == "factor") {
    x = as.character(x)
    classx = "character"
  }
  if (classx == "character") {
    x[x==""] = NA
  }
  i0 = is.na(x)
  nmiss = sum(i0)
  iv = which(!i0)
  xv = x[iv]
  return(xv)
}

pred_profile = function(x) {
  n = length(x)
  xv = preprocess(x)
  classx = class(xv)[1]
  nmiss = n - length(xv)
  nd = length(unique(xv)) + (nmiss > 0)
  if (classx == "numeric") {
    if (nmiss < n) {
      if (quantile(xv - floor(xv), 0.99) < 1e-3) {cl = "int"}
      else {cl = "cont"}
    } else {cl = "int"}
  } else if (classx == "integer") { 
    cl = "int" 
  } else if (classx == "logical") {
    cl = "int"
  } else {
    cl = classx[1]
  }
  if (cl == "cont" || cl == "int") {
    qx = c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1)
    qq = quantile(xv, qx)
    pprofile = list(
      flag_predictor = ifelse(nd > 1, 1, 0),
      cl = cl,
      n = n,
      nmiss = nmiss,
      pmiss = nmiss / n,
      nd = nd,
      quantiles = tibble(p = qx, x = qq))
  } else if (cl == "character") {
    ut = table(xv)
    if (nmiss>0) {ut["NA"] = nmiss}
    pprofile = list(
      flag_predictor = ifelse(nd > 1, 1, 0),
      cl = cl,
      n = n,
      nmiss = nmiss,
      pmiss = nmiss / n,
      nd = nd,
      quantiles = tibble(p = cumsum(as.integer(ut)), x = names(ut)))
  } else {
    pprofile = list(
      flag_predictor = 0,
      cl = cl,
      n = n,
      nmiss = nmiss,
      pmiss = nmiss / n,
      nd = nd + ifelse(nmiss > 0, 1, 0))
  }
  return(pprofile)
}

grid_decode = function(code, cl) {
  if (cl %in% c("cont","int")) {
    tokens = strsplit(code, split = ", ")[[1]]
    if (length(tokens) == 1 && tokens == "NA") {
      borders = NA
    } else {
      borders = as.numeric(tokens)
    }
  } else if (cl == "character") {
    l1tokens = strsplit(code, split = "\", \"")
    if (length(l1tokens)>0) {
      l2df = list()
      for (i in 1:length(l1tokens[[1]])) {
        l = l1tokens[[1]][i]
        l2tokens = strsplit(l, "\" \"")
        if (length(l2tokens)>0) {
          ll = gsub('"', '', l2tokens[[1]]) # remove quotation marks
          l2df[[i]] = data.frame(level = ll, bin = i, stringsAsFactors = FALSE)
        }
      }
      borders = bind_rows(l2df)
    } else {
      borders = data.frame(level = "", bin = 1)
    }
  } else {
    borders = NULL
  }
  return(borders)
}

grid_encode = function(borders, cl) {
  if (cl %in% c("cont","int")) {
    if (all(!is.na(borders))) {
      code = paste(borders, collapse = ", ")
    } else {
      code = "NA"
    }
  } else if (cl == "character") {
    if (nrow(borders)>=1) {
      code = paste0('"', borders$level[1], '"')
      if (nrow(borders)>=2) {
        for (i in 2:nrow(borders)) {
          if (borders$bin[i] != borders$bin[i-1])
            code = paste0(code, ", \"", borders$level[i], '"')
          else
            code = paste0(code, " \"", borders$level[i], '"')
        }
      }
    } else {
      code = ""
    }
  } else {
    code = ""
  }
  return(code)
}

make_summary = function(d, t0) {
  tdata = d %>% group_by(id_bin) %>% 
    summarise(
      n = length(id_bin),
      m = sum(target, na.rm=TRUE),
      r = round(100 * length(id_bin) / nrow(d), 2),
      p = round(100*mean(target, na.rm=TRUE), 2))
  tdata$WoE = woeagg(tdata)
  t = t0 %>% left_join(tdata, by = "id_bin") %>% select(id_bin, bin, n, m, r, p, WoE)
  # assign correct values to empty bins
  t$n[is.na(t$n)] = 0
  t$m[t$n == 0] = 0
  t$r[t$n == 0] = 0
  t$p[t$n == 0] = round(100*mean(d$target, na.rm=TRUE), 2)
  t$WoE[t$n == 0] = 0
  return(t)  
}

binning_manual = function(x, target, itest = NULL, pprofile, code, missbin) {
  i0 = is.na(x)
  xv = preprocess(x)
  if (code == "" || code == "NA" || pprofile$flag_predictor == 0) { # not eligible for binning
    b = pprofile
    b$flag_predictor = 0
    b$nbins = 1
    b$missbin = 1
    b$borders = NA
    b$code = code
    b$gini_test = 0
    b$lift_test = 1
    b$gini_train = 0
    b$lift_train = 1
    return(b)
  }
  if (pprofile$cl %in% c("cont", "int")) {
    h0 = grid_decode(code, pprofile$cl)
    h = c(-Inf, unique(h0), +Inf)
    yv = cut(as.numeric(xv), breaks = h, right = FALSE)
    nbin0 = length(levels(yv))
    if (missbin <= nbin0) { # missing is merged with an existing bin
      levels(yv)[missbin] = paste0(levels(yv)[missbin], ", NA")
    } else { # missing makes a standalone bin K+1
      missbin = nbin0 + 1
      levels(yv)[missbin] = "NA"
    }
    grp = rep(0,pprofile$n)
    grp[!i0] = as.integer(yv)
    grp[i0] = missbin
    nbins = length(levels(yv))
    t0 = tibble(bin = levels(yv), id_bin = union(1:nbins, missbin))
    d = tibble(id_bin = grp, target = target)
    if (is.null(itest)) {
      sn = make_summary(d, t0)
      gini_train = gini_table(sn)
      lift_train = max(sn$p/100) / mean(target)
      st = NULL
      gini_test = NA
      lift_test = NA
    } else {
      sn = make_summary(d[-itest,], t0)
      gini_train = gini_table(sn)
      lift_train = max(sn$p/100) / mean(target[-itest])
      st = make_summary(d[itest,], t0)
      st$WoE = sn$WoE # WoE musi zustat zafixovane v train
      gini_test = gini_table(st)
      lift_test = max(st$p/100) / mean(target[itest])
    }
    b = pprofile 
    b$nbins = nbins
    b$missbin = missbin
    b$borders = h0
    b$code = code
    b$summary_train = sn
    b$summary_test = st
    b$gini_test = gini_test
    b$lift_test = lift_test
    b$gini_train = gini_train
    b$lift_train = lift_train
  } else if (pprofile$cl == "character") {
    # prevodni tabulka level -> bin
    tt = grid_decode(code, pprofile$cl)
    xb = rep(0, pprofile$n)
    nbins = nrow(tt)
    grp = rep(0, pprofile$n)
    for (j in 1:nbins) {
      ij = which(x == tt$level[j]) #TODO nebude fungovat kdyz se merguji biny
      xb[ij] = tt$bin[j]
      grp[ij] = j
    }
    xb[i0] = missbin
    t0 = tibble(id_bin = 1:nbins, bin = tt$level) %>% arrange(bin)
    # bin labels - ted neni funkcni, bny se nemerguji
    # t0$bin = ""
    # for (k in 1:nrow(t0)) {
    #   t0$bin[k] = paste(unique(tt$level[tt$bin==k]), collapse = " ")
    # }
    d = tibble(id_bin = grp, target = target)
    if (is.null(itest)) {
      sn = make_summary(d, t0)
      gini_train = gini_table(sn)
      lift_train = max(sn$p/100) / mean(target)
      st = NULL
      gini_test = NA
      lift_test = NA
    } else {
      sn = make_summary(d[-itest,], t0)
      gini_train = gini_table(sn)
      lift_train = max(sn$p/100) / mean(target[-itest])
      st = make_summary(d[itest,], t0)
      st$WoE = sn$WoE # WoE musi zustat zafixovane v train
      gini_test = gini_table(st)
      lift_test = max(st$p/100) / mean(target[itest])
    }
    b = pprofile 
    b$nbins = nbins
    b$missbin = missbin
    b$borders = sn$bin
    b$code = code
    b$summary_train = sn
    b$summary_test = st
    b$gini_test = gini_test
    b$lift_test = lift_test
    b$gini_train = gini_train
    b$lift_train = lift_train
  } else {
    b = pprofile
  }
  return(b)
}

woeagg = function(sumt) {
  # rows of sumt are levels of categorical variable
  # sumt$n  absolute number of observations for each level
  # sumt$m  absolute number of events for each level
  n1 = sumt$m
  n0 = sumt$n - sumt$m
  sigj = log(n1/n0)
  sig = log(sum(n1)/sum(n0))
  WoE = sigj - sig
  WoE[n1==0] = -3
  WoE[n0==0] = +3
  return(round(WoE, digits = 3))
}

gini_table = function(tt) {
  if (nrow(tt)>1) {
    # columns n, m, p needed
    t2 = tt[order(tt$WoE),]
    t2$F0 = cumsum(t2$n - t2$m)/sum(t2$n - t2$m)
    t2$F1 = cumsum(t2$m)/sum(t2$m)
    t2$gi = 0 # allocate new column
    for (j in 2:nrow(t2)) {
      t2$gi[j] = t2$F1[j]*t2$F0[j-1] - t2$F0[j]*t2$F1[j-1]
    }
    sum(t2$gi)
  } else {
    return(0)
  }
}

suggest_decimal = function(pprofile) {
  if (!is.null(pprofile$quantiles)) {
    dec = 2
  } else if (pprofile$cl=="int") { 
    dec = 0 
  } else if (pprofile$cl=="cont") {
    p25 = pprofile$quantiles %>% filter(p = 0.25) %>% select(x)
    p75 = pprofile$quantiles %>% filter(p = 0.75) %>% select(x)
    if (p75 > p25) {
      dec = round(3 - log10(p75 - p25), 0)
    } else {
      dec = 2
    }
  }
  unname(dec)
}

bins_eq_frequency = function(x, nbins, dec) {
  pb = (1:(nbins-1))/nbins
  h0 = quantile(x, pb, na.rm = TRUE, names = FALSE, type = 3)
  h1 = unique(round(h0, dec))
  if (sum(x < h1[1], na.rm = TRUE) == 0) { # mass probability at min value greater than 1/nbins
    if (length(h1) >= 2) {
      h2 = h1[2:length(h1)]
    } else if (max(x, na.rm = TRUE) == h1) { # mass probability at min value equal to 1 (all observations concentrated at one value)
      h2 = NA
    } else { # mass probability at min value greater than (nbins - 1)/nbins
      h2 = min(x[x > h1])
    }
  } else {
    h2 = h1
  }
  return(h2)
}

binning_eq_frequency = function(x, target, itest = NULL, pprofile, nbins) {
  xv = preprocess(x)
  if (pprofile$flag_predictor == 1 && pprofile$cl %in% c("cont","int")) {
    dec = suggest_decimal(pprofile)
    borders = bins_eq_frequency(xv, nbins, dec)
    code = grid_encode(borders, pprofile$cl)
    missbin = length(borders) + 2
  } else if (pprofile$flag_predictor == 1 && pprofile$cl == "character") {
    nbins = pprofile$nd
    v = sort(unique(xv))
    borders = tibble(level = v, bin = 1:(nbins - (pprofile$nmiss > 0)), stringsAsFactors = FALSE)
    code = grid_encode(borders, pprofile$cl)
    missbin = nbins
  } else {
    code = ""
    missbin = 1
  }
  b = binning_manual(x, target, itest, pprofile, code, missbin)
  return(b)
}

bins_ctree = function(x, target, p = 0.05) {
  df = tibble(x = x, target = target)
  ctree = partykit::ctree(
    target ~ x,
    data = df,
    na.action = na.exclude,
    control = partykit::ctree_control(minbucket = ceiling(round(p*nrow(df)))))
  bins = partykit::width(ctree)
  if (bins<2) {return(NA)}
  # Append cutpoinstop()ts in a table (Automated)
  cutvct = data.frame(matrix(ncol=0,nrow=0)) # Shell
  n = length(ctree) # Number of nodes
  for (i in 1:n) {
    cutvct = rbind(cutvct, ctree[i]$node$split$breaks)
  }
  cutvct = cutvct[order(cutvct[,1]),] # Sort / converts to a ordered vector (asc)
  cutvct = ifelse(cutvct<0,trunc(10000*cutvct)/10000,ceiling(10000*cutvct)/10000) # Round to 4 dec. to avoid borderline cases
  return(cutvct)
}

binning_ctree = function(x, target, itest = NULL, pprofile) {
  i0 = is.na(x)
  xv = preprocess(x)
  if (pprofile$flag_predictor == 1 && pprofile$cl %in% c("cont","int")) {
    borders = bins_ctree(xv, target[!i0])
    #dec = suggest_decimal(pprofile)
    #borders = unique(round(borders, dec))
    code = grid_encode(borders, pprofile$cl)
    missbin = length(borders) + 2
  } else if (pprofile$flag_predictor == 1 && pprofile$cl == "character") {
    nbins = pprofile$nd
    v = sort(unique(xv))
    borders = tibble(level = v, bin = 1:nbins, stringsAsFactors = FALSE)
    code = grid_encode(borders, pprofile$cl)
    missbin = nbins
  } else {
    code = ""
    missbin = 1
  }
  b = binning_manual(x, target, itest, pprofile, code, missbin)
  return(b)
}

pred_binning = function(x, target, itest = NULL, pprofile, method = "EqualFrequency", ...) {
  if (method == "EqualFrequency") {
    binning_eq_frequency(x, target, itest, pprofile, ...)
  } else if (method == "ctree") {
    binning_ctree(x, target, itest, pprofile, ...)
  } else {
    stop(sprintf("Binning method %s is currently not supported.\n", method))
  }
}

autobinning = function(tbl, target, predictors, itest = NULL, method = "EqualFrequency", verbose = TRUE, maxlength = 32, ...) {
  # itest - indicies of testing set
  if (length(target) == 1) { # target is either column index or column name
    tg = tbl[[target]]
  } else {
    tg = target
  }
  tn = names(tbl)
  if (class(predictors) == "character") {
    ipred = which(tn %in% predictors)
  } else {
    ipred = predictors
  }
  np = length(ipred)
  woenames = woe_names(tn[ipred], maxlength = maxlength, verbose = TRUE)
  L = tibble(
    ordnum = 1:np, 
    varnum = ipred, 
    varname = tn[ipred], 
    class = "",
    flag_predictor = 0,
    pmiss = 0,
    nd = 0,
    nbins = 0,
    gini_test = 0,
    lift_test = 0,
    gini_train = 0,
    lift_train = 0,
    missbin = 0,
    code = "",
    woename = woenames,
    binerr = 0,
    quantiles = vector("list", np),
    binning_train = vector("list", np),
    binning_test = vector("list", np),
    PSI = 0)
  if (is.null(itest)) {
    dummy_train = tibble(id_bin = 1, bin = "all", n = length(tg), m = sum(tg), r = 100, p = 100*mean(tg), WoE = 0)
    dummy_test = tibble(id_bin = 1, bin = "test", n = 0, m = 0, r = NA, p = NA, WoE = NA)
  } else {
    dummy_train = tibble(id_bin = 1, bin = "train", n = length(tg) - length(itest), m = sum(tg[-itest]), r = 100, p = 100*mean(tg[-itest]), WoE = 0)
    dummy_test = tibble(id_bin = 1, bin = "test", n = length(itest), m = sum(tg[itest]), r = 100, p = 100*mean(tg[itest]), WoE = 0)
  }
  for (i in 1:np) {
    j = ipred[i]
    if (verbose == TRUE) {
      cat(sprintf("%04d/%04d #%04d %-32s ", i, np, j, str_trunc(tn[j], 32)))
    }
    xi = tbl[[j]]
    pp = pred_profile(xi)
    err_msg = ""
    tryCatch({ 
      bb = pred_binning(xi, tg, itest, pp, method, ...)
      bb$name = tn[j]
      L$class[i] = bb$cl
      L$flag_predictor[i] = bb$flag_predictor
      L$pmiss[i] = bb$pmiss
      L$nd[i] = bb$nd
      L$nbins[i] = bb$nbins
      L$gini_test[i] = bb$gini_test
      L$lift_test[i] = bb$lift_test
      L$gini_train[i] = bb$gini_train
      L$lift_train[i] = bb$lift_train
      L$missbin[i] = bb$missbin
      L$code[i] = bb$code
      L$quantiles[[i]] = bb$quantiles
      if (!is.null(bb$summary_train)) {L$binning_train[[i]] = bb$summary_train}
      else {L$binning_train[[i]] = dummy_train}
      if (!is.null(bb$summary_test)) {L$binning_test[[i]] = bb$summary_test}
      else {L$binning_test[[i]] = dummy_test}
      L$PSI[i] = PSI(L$binning_train[[i]], L$binning_test[[i]]) # population stability index
      if (verbose == TRUE) {
        cat(sprintf(" cl = %-10s, pmiss = %.2f, nd = %5d, nbins = %3d, gini = %.3f, lift = %.2f\n", 
          bb$cl, bb$pmiss, bb$nd, bb$nbins, bb$gini_test, bb$lift_test))
      }
    }, error = function(e) { # <<- assignment is necessary to evaluate in parent namespace
      L$binerr[i] = 1
      #L$quantiles[[i]] <<- pp$quantiles
      if (verbose == TRUE) {
        cat(sprintf(" cl = %-10s, pmiss = %.2f, nd = %5d, error=%s\n", 
          pp$cl, pp$pmiss, pp$nd, e))
      }
    })
  }
  return(L)
}

pred_grp = function(x, b) {
  i0 = is.na(x)
  xv = preprocess(x)
  if (b$code == "" || b$flag_predictor == 0) { # not eligible for binning
    grp = rep(0, length(x))
    return(grp)
  }
  if (b$class %in% c("cont", "int")) {
    h0 = grid_decode(b$code, b$class)
    h = c(-Inf, unique(h0), +Inf)
    y = cut(as.numeric(xv), breaks = h, right = FALSE)
    nbin0 = length(levels(y))
    if (b$missbin <= nbin0) { # missing is merged with an existing bin
      levels(y)[b$missbin] = paste0(levels(y)[b$missbin], ", NA")
      y[i0] = levels(y)[b$missbin]
      grp = as.integer(y)
    } else { # missing makes a standalone bin K+1
      missbin = nbin0 + 1
      levels(y)[missbin] = "NA"
      w=y
      w[i0] = "NA"
      w[!i0] = y
      grp = as.integer(w)
    }
  } else if (b$class == "character") {
    # tt je prevodni tabulka level -> bin
    tt = grid_decode(b$code, b$class)
    grp = rep(0, length(x))
    for (j in 1:nrow(tt)) {
      ij = which(x == tt$level[j])
      grp[ij] = j
    }
    grp[i0] = b$missbin
  } else {
    grp = rep(0, length(x))
  }
  return(grp)
}

pred_woe = function(x, b) {
  gg = pred_grp(x, b)
  woe = rep(0, length(x))
  if (!is.null(b$binning_train)) {
    tt = b$binning_train[[1]]
    for (j in 1:nrow(tt)) {
      ij = which(gg == j)
      woe[ij] = tt$WoE[j]
    }
  }
  return(woe)
}

tbl_woe = function(tbl, tbl_binning, keep = c(), verbose = TRUE) {
  np = nrow(tbl_binning)
  n = nrow(tbl)
  bnm = tbl_binning$varname
  nm = union(intersect(keep, names(tbl)), bnm)
  woe = tbl[nm]
  for (i in 1:np) {
    nmi = tbl_binning$varname[i]
    j = which(nm == nmi)
    if (verbose == TRUE) {
      cat(sprintf("WoE transformation %04d/%04d #%04d %-32s\n", i, np, j, nmi))
    }
    b = tbl_binning[i,]
    woe[[j]] = pred_woe(tbl[[nmi]], b)
    names(woe)[j] = tbl_binning$woename[i]
  }
  return(woe)
}

write_sql = function(outfile, intable, outtable, tbl_binning, keep = c()) {
# Compiles SQL file which calculates WoE transformation of predictors from intable based on array of binning objects (tbl_binning).
# outfile       path to the output SQL file
# intable       name of the input table to be used in FROM clause of the SQL statement
# outtable      name of the output table to be use in INTO clause of the SQL statement
# tbl_binning   array of binning objects. Each binning object b contains
#   $name              code name of the predictor (without woe_ prefix)
#   $flag_predictor    boolean value if the field qualifies as a predictor
#   $cl                class, one of character, int, cont
#   $summary           data frame with bins and their WoEs
#   $borders           vector of numerical bin borders
# keep          names of the columns to be copied from input table without WoE transformation (ids, observation date, target, ...)
  fid = file(outfile, "w")
  cat("-- Weight of evidence SQL generated by autobinning\n", file = fid)
  cat(sprintf("-- creates table %s from table %s\n", outtable, intable), file = fid)
  cat(sprintf("if object_id('%s') is not null drop table %s\ngo\n\n", outtable, outtable), file = fid)
  cat("select ", file = fid)
  prefix = ""
  if (length(keep)>0) {
    for (i in 1:length(keep)) {
      cat(sprintf("%s%s", prefix, keep[i]), file = fid)
      prefix = "\n  , "
    }
  }
  for (i in 1:length(tbl_binning)) {
    b = tbl_binning[[i]]
    if (b$flag_predictor == 1) {
      cat(sprintf("%scase", prefix), file = fid)
      tt = b$summary
      for (j in 1:nrow(tt)) {
        if (tt$bin[j] == "NA") {
          cat(sprintf("\n    when %s is null then %.03f", b$name, tt$WoE[j]), file = fid)
        } else if (b$cl == "character") {
          cat(sprintf("\n    when %s = '%s' then %.03f", b$name, tt$bin[j], tt$WoE[j]), file = fid)
        } else if (b$cl %in% c("cont", "int")) {
          if (j == 1) {
            cat(sprintf("\n    when %s < %g then %.03f -- interval %s", b$name, b$borders[1], tt$WoE[j], tt$bin[j]), file = fid)
          } else if (j < nrow(tt) && tt$bin[j+1] != "NA") {
            cat(sprintf("\n    when %s >= %g and %s < %g then %.03f -- interval %s", b$name, b$borders[j-1], b$name, b$borders[j], tt$WoE[j], tt$bin[j]), file = fid)
          } else {
            cat(sprintf("\n    when %s >= %g then %.03f -- interval %s", b$name, b$borders[j-1], tt$WoE[j], tt$bin[j]), file = fid)
          }
        }
      }
      cat(sprintf("\n    else 0.000 end woe_%s", b$name), file = fid)
      prefix = "\n  , "
    }
  }
  cat(sprintf("\ninto %s\nfrom %s\ngo\n", outtable, intable), file = fid)
  close(fid)
}

woe_names = function(name, maxlength = 32, verbose = TRUE) {
# modifies variable names to woe_* to a specified number of characters
  wname = rep("", length(name))
  stemmax = maxlength - 4 # without woe_
  for (i in 1:length(name)) {
    if (nchar(name[i]) <= stemmax) { # woe_XXXX does not exceed 32 characters
      wname[i] = sprintf("woe_%s", name[i])
    } else {
      truncname = sprintf("woe_%s", substr(name[i], 1, stemmax))
      if (truncname %in% wname) { # we will add numbers at the end
        k = 1
        while (sprintf("woe_%s%d", substr(name[i], 1, stemmax - 1 - trunc(log10(k))), k) %in% wname) {
          k = k + 1
        }
        wname[i] = sprintf("woe_%s%d", substr(name[i], 1, stemmax - 1 - trunc(log10(k))), k)
      } else {
        wname[i] = truncname
      }
      if (verbose == TRUE) {
        cat(sprintf("name #%d %s truncated to %s\n", i, name[i], wname[i]))
      }
    }
  }
  return(wname)
}

write_sas = function(outfile, intable, outtable, tbl_binning, keep = c(), verbose = TRUE) {
# Compiles SAS file with data step which calculates WoE transformation of predictors from intable based on array of binning objects (tbl_binning).
# outfile       path to the output SAS file
# intable       name of the input table to be used in SET clause of the data step
# outtable      name of the output table to be use in DATA clause of the data step
# tbl_binning   array of binning objects. Each binning object b contains
#   $name              code name of the predictor (without woe_ prefix)
#   $flag_predictor    boolean value if the field qualifies as a predictor
#   $cl                class, one of character, int, cont
#   $summary           data frame with bins and their WoEs
#   $borders           vector of numerical bin borders
# keep          names of the columns to be copied from input table without WoE transformation (ids, observation date, target, ...)
# verbose       boolean indicator of woe_ variable name truncations should be reported
  fid = file(outfile, "w")
  cat("/* Weight of evidence SAS datastep generated by autobinning */\n", file = fid)
  cat(sprintf("/* creates table %s from table %s */\n\n", outtable, intable), file = fid)
  cat(sprintf("data %s (keep = woe_:", outtable), file = fid)
  if (length(keep)>0) {
    for (i in 1:length(keep)) {
      cat(sprintf(" %s", keep[i]), file = fid)
    }
  }
  cat(sprintf(");\nset %s;\n", intable), file = fid)
  nm = sapply(tbl_binning, function(x) {getElement(x, "name")})
  wname = woe_names(nm, maxlength = 32, verbose = verbose)
  for (i in 1:length(tbl_binning)) {
    b = tbl_binning[[i]]
    if (b$flag_predictor == 1) {
      cat(sprintf("\nselect; /* %s */", b$name), file = fid)
      tt = b$summary
      for (j in 1:nrow(tt)) {
        if (tt$bin[j] == "NA") {
          cat(sprintf("\n  when (missing(%s)) %s = %.03f;", b$name, wname[i], tt$WoE[j]), file = fid)
        } else if (b$cl == "character") {
          cat(sprintf("\n  when (%s = '%s') %s = %.03f;", b$name, tt$bin[j], wname[i], tt$WoE[j]), file = fid)
        } else if (b$cl %in% c("cont", "int")) {
          if (j == 1) {
            cat(sprintf("\n  when (%s < %g) %s = %.03f; /* interval %s */", b$name, b$borders[1], wname[i], tt$WoE[j], tt$bin[j]), file = fid)
          } else if (j < nrow(tt) && tt$bin[j+1] != "NA") {
            cat(sprintf("\n  when (%g <= %s < %g) %s = %.03f; /* interval %s */", b$borders[j-1], b$name, b$borders[j], wname[i], tt$WoE[j], tt$bin[j]), file = fid)
          } else {
            cat(sprintf("\n  when (%s >= %g) %s = %.03f; /* interval %s */", b$name, b$borders[j-1], wname[i], tt$WoE[j], tt$bin[j]), file = fid)
          }
        }
      }
      cat(sprintf("\n  otherwise %s = 0.000;\nend;", wname[i]), file = fid)
    }
  }
  cat(sprintf("\nrun;\n", outtable, intable), file = fid)
  close(fid)
}

compatible_axes = function(y1, y2) {
  # suggests two vertical axes for plotting ranges [0,y1] and [0,y2]
  # The ylim is selected from the sequence 0.5,1,2.5,5,10,25,50 such that
  # y1 is between (40%, 100%] of ylim1 and y2 is between (40%, 100%] of ylim2.
  # each axis can then be split into five equal bins and labeled with six tick marks,
  # e.g. 0, 5, 10, 15, 20, 25.
  c = c(1, 2.5, 5, 10) # one cycle
  lc = log10(c)
  y = c(y1, y2)
  ar = log10(y) # aspect ratio
  base = floor(ar)
  rem = ar - base # remainder
  ia = c(0, 0)
  for (i in 1:2) {
    if (rem[i]==lc[1]) ia[i] = 1
    else if (rem[i]<lc[2]) ia[i] = 2
    else if (rem[i]<lc[3]) ia[i] = 3
    else ia[i] = 4
  }
  ax = 10^base * c[ia]
}

plotbarline = function(ybar, yline, xlabels, xname, title) {
  color_ybar = rgb(0/255,151/255,169/255)
  color_yline = rgb(237/255,139/255,0/255)
  Gray50 = rgb(0.5,0.5,0.5)
  Gray70 = rgb(0.7,0.7,0.7)
  n = length(ybar)
  if (length(yline)!=n) stop("ybar and yline must be of the same length")
  if (length(xlabels)!=n) stop("xlabels must be of the same length as ybar and yline")
  ax = compatible_axes(max(ybar), max(yline))
  xa = 1:n
  ya0 = c(0,0.2,0.4,0.6,0.8,1)
  ya1 = ax[1]*ya0
  ya2 = ax[2]*ya0
  yl1 = as.character(ya1)
  yl2 = as.character(ya2)
  barlabels = as.character(ybar)
  linelabels = as.character(yline)
  plot.new()
  plot.window(xlim=c(0.5,n+0.5), ylim=c(0,1))
  # bottom axis
  ba = axis(side=1, at = xa, labels = xlabels, tick = TRUE, line = NA,
            pos = NA, outer = FALSE, font = NA, lty = "solid",
            lwd = 1, lwd.ticks = 1, col = Gray50, col.ticks = Gray50,
            hadj = NA, padj = NA)
  # left axis
  la = axis(side=2, at = ya0, labels = yl1, tick = TRUE, line = NA,
            pos = NA, outer = FALSE, font = NA, lty = "solid",
            lwd = 1, lwd.ticks = 1, col = Gray50, col.ticks = Gray50,
            hadj = NA, padj = NA)
  # right axis
  ra = axis(side=4, at = ya0, labels = yl2, tick = TRUE, line = NA,
            pos = NA, outer = FALSE, font = NA, lty = "solid",
            lwd = 1, lwd.ticks = 1, col = Gray50, col.ticks = Gray50,
            hadj = NA, padj = NA)
  abline(h=ya0, col=Gray70)
  ybaradj = ybar / ax[1]
  ylineadj = yline / ax[2]
  rect(xleft=xa-0.4, ybottom=rep(0,n), xright=xa+0.4, ytop=ybaradj,
       col = color_ybar, border = color_ybar, lty = "solid", lwd = 1)
  text(xa, ybaradj, barlabels, cex=1, pos=3, col=color_ybar)
  lines(x=xa, y=ylineadj, col = color_yline, type="l", lwd=3)
  lines(x=xa, y=ylineadj, col = color_yline, type="p", lwd=1, pch=15, cex=1.5)
  text(xa, ylineadj, linelabels, cex=1, pos=3, col=color_yline)
  title(main=title, xlab=xname)
  #legend(1, 0.2, legend = c("event rate","count"))
}

plot_binning = function(binning, path = ".") {
  for (i in 1:nrow(binning)) {
    binning_row = binning[i,]
    fname = sprintf("%s/%s.png", path, binning_row$varname)
    cat(sprintf("#%04d %s\n", i, fname))
    png(fname, width = 600, height = 400)
    plot_binning_train(binning_row)
    dev.off()
  }
}

plot_binning_train = function(binning_row) {
  btab = binning_row$binning_train[[1]]
  plotbarline(
    ybar = btab$n, 
    yline = btab$p, 
    xlabels = btab$bin, 
    xname = binning_row$varname,
    title = sprintf("%.3f on train (%d observations)", binning_row$gini_train, sum(btab$n)))
}

plot_binning_test = function(binning_row) {
  btab = binning_row$binning_test[[1]]
  plotbarline(
    ybar = btab$n, 
    yline = btab$p, 
    xlabels = btab$bin, 
    xname = binning_row$varname,
    title = sprintf("%.3f on test (%d observations)", binning_row$gini_test, sum(btab$n)))
}

adjust_binning = function(L, tbl, itest = NULL, ni, target, code, missbin) {
  # adjust_binning pouzivej po autobinningu
  # priklad:
  # b = adjust_binning(
  #   b, #cely autobinning
  #   dt, #input data
  #   its, # test indicies
  #   "team_churn_dif_12M", #jeden prediktor zmenim
  #   target = "target", 
  #   code = "0, 0.07, 0.11, 0.18, 0.32",
  #   missbin = 7)
  i = which(L$varname == ni)
  yi = L[i,]
  xi = tbl[[ni]]
  pp = pred_profile(xi)
  bi = binning_manual(xi, tbl[[target]], itest, pp, code = code, missbin)
  yi$class = bi$cl
  yi$flag_predictor = bi$flag_predictor
  yi$pmiss = bi$pmiss
  yi$nd = bi$nd
  yi$nbins = bi$nbins
  yi$gini_train = bi$gini_train
  yi$gini_test = bi$gini_test
  yi$lift_train = bi$lift_train
  yi$lift_test = bi$lift_test
  yi$missbin = bi$missbin
  yi$code = bi$code
  yi$binning_train[[1]] = bi$summary_train
  yi$binning_test[[1]] = bi$summary_test
  L[i,] = yi
  return(L)
}

PSI = function(tt1, tt2) {
  tt = left_join(
    tt1 %>% filter(n>0) %>% rename(n1 = n), 
    tt2 %>% rename(n2 = n), 
    by = "id_bin")
  p_expect = tt$n1/sum(tt$n1)
  p_actual = tt$n2/sum(tt$n2, na.rm = TRUE)
  izpa = which(p_actual == 0)
  if (length(izpa)>0) {
    lae = log(p_actual[-izpa] / p_expect[-izpa])
    PSI = (p_actual[-izpa] - p_expect[-izpa]) * lae + 0.25*length(izpa)
  } else {
    lae = log(p_actual / p_expect)
    PSI = (p_actual - p_expect) * lae
  }
  return(sum(PSI))
}

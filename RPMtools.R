# Validation tools
knitr::opts_chunk$set(
  warning = F, # show warnings
  message = F, # show messages
  echo = F  # show R code
)
library(tidyverse)
library(odbc)
library(dbplyr)
library(data.table)
library(CatEncoders)
library(glue)
library(plotly)
library(htmltools)
library(dplyr)
library(sf)
library(gridExtra)
library(tidyr)
library(lubridate)
library(reshape2)
library(reticulate)
library(ggplot2)
library(ParBayesianOptimization)
library(mlbench)
library(recipes)
library(resample)
library(xgboost)
library(caret)
library(Matrix)
library(magrittr)

# library(SHAPforxgboost)
# con <- dbConnect(odbc::odbc(), "AUTOSQL02" ,database = "MISAnalystDB")
# con <- dbConnect(odbc::odbc(), "AUTOSQL02" ,database = "MISReportsDB")
# snapshot1 <- con  %>% tbl(in_schema("dbo" , "GenericBackfill_2018Partial"))
# snapshot2 <- con  %>% tbl(in_schema("dbo" , "GenericBack"))


options(scipen=999)


KT_create_equal_bin = function(weight,nbin){
  cumulative_sum = cumsum(weight)
  bins =  cut(cumulative_sum, breaks = nbin , labels = F)
  return(bins)
}

######################## Gini ############################

# https://www.kaggle.com/c/liberty-mutual-fire-peril/discussion/9880 Kaggle competition modelling fire loss cost

KT_calc_gini_old = function(actual, weight, predicted){
  df = data.frame(actual = actual, weight = weight, predicted =
                    predicted)
  # create random number sort key so ties will be resolved in random order
  k = length(df$actual)
  df$rkey = runif(k)
  df = df[order(df$predicted, df$rkey),]
  df$random = cumsum((df$weight/sum(df$weight)))
  totalPositive = sum(df$actual * df$weight)
  df$cumPosFound = cumsum(df$actual * df$weight)
  df$Lorentz = df$cumPosFound / totalPositive
  n = nrow(df)
  gini = sum(df$Lorentz[-1]*df$random[-n]) - sum(df$Lorentz[-n]
                                                 * df$random[-1])
  return(gini)
}


KT_calc_gini <- function(actual, weight,predicted){
  df= data.frame(actual,predicted,weight)
  actual = as.numeric(df$actual)
  weight = as.numeric(df$weight)
  s_idx = order(df$predicted)
  w_s = weight[s_idx]
  a_s = actual[s_idx]
  a_c = cumsum(a_s*w_s)
  w_c = cumsum(w_s)
  a_c = a_c/a_c[length(a_c)]
  w_c = w_c/w_c[length(w_c)]
  return(1-2*pracma::trapz(w_c,a_c) )}

KT_calc_gini_norm =function(actual, weight, predicted){
  return(KT_calc_gini(actual, weight, predicted)/KT_calc_gini(actual, weight, actual))
}


KT_resample_gini = function(n,actual, weight, predicted, normalize = FALSE){
  gini_vector = numeric()
  df = data.frame(actual = actual, weight = weight, predicted =predicted)
  for (x in sample(1:n*33, n, replace=FALSE) ){
    set.seed(x)
    test = df[sample(nrow(df), size = nrow(df), replace = T), ]
    if (normalize ==F){
      gini_vector = c(gini_vector , KT_calc_gini( actual =  test$actual, weight =  test$weight , predicted =  test$predicted))
    }
    else {
      gini_vector= c(gini_vector , KT_calc_gini_norm( actual =  test$actual, weight =  test$weight , predicted =  test$predicted))
    }
  }
  return(gini_vector)
}

KT_plot_compare_gini =  function(n,actual, weight, base, challenger, normalize = FALSE){
  base_gini = KT_resample_gini(n,actual, weight, base, normalize = normalize)
  challenger_gini = KT_resample_gini(n,actual, weight, challenger, normalize = normalize) 
  gini_df = cbind( data.frame( challenger_gini =  challenger_gini),
                   data.frame(   base_gini = base_gini)) 
  gini_df %>% mutate(challenger_win_ind = ifelse(challenger_gini>base_gini , 1,0))  -> gini_df
  challenger_win_rate = mean(gini_df$challenger_win_ind)
  gini_df %>% select(-challenger_win_ind) %>%
    melt() %>%
    ggplot(.,aes(x = value, colour = variable , fill = variable)  ) + geom_density(alpha = 0.3) +
    ggtitle(glue("Gini comparison | Challenger win rate = {challenger_win_rate}"))
  
}


KT_plot_compare_3models_gini =  function(n,actual, weight, base, challenger,challenger2, normalize = FALSE){
  base_gini = KT_resample_gini(n,actual, weight, base, normalize = normalize)
  challenger_gini = KT_resample_gini(n,actual, weight, challenger, normalize = normalize) 
  challenger2_gini = KT_resample_gini(n,actual, weight, challenger2, normalize = normalize) 
  gini_df = cbind( 
                   data.frame(   base_gini = base_gini),
                   data.frame( challenger_gini =  challenger_gini),
                   data.frame( challenger2_gini =  challenger2_gini)) 
  gini_df %>% mutate(challenger_win_ind = ifelse(challenger_gini>base_gini , 1,0))  -> gini_df
  challenger_win_rate = mean(gini_df$challenger_win_ind)
  gini_df %>% select(-challenger_win_ind) %>%
    melt() %>%
    ggplot(.,aes(x = value, colour = variable , fill = variable)  ) + geom_density(alpha = 0.3) +
    ggtitle(glue("Gini comparison"))
  
}

######################## R^2 ############################

KT_Rsq = function(actual,pred){
  r_sq = 1- (sum((actual-pred)^2)/sum((actual-mean(actual))^2))
  return( r_sq)
}
KT_weighted_Rsq = function(actual,pred, weight){
  r_sq = 1- (sum(((actual-pred)^2)*weight)/sum(((actual-mean(actual))^2)*weight))
  return( r_sq)
}



######################## Lift ############################

KT_calc_lift = function(pred , actual, weight ,nbin){
  pred =   pred*(sum(actual)/sum(pred*weight )) # rebase
  lift_df = data.frame(pred , actual, weight )
  lift_df = lift_df %>% 
    filter(weight>0) %>%
    arrange(pred) %>%
    mutate(pred=pred*weight,
           bin = KT_create_equal_bin(weight,nbin) )
  lift_df_agg= lift_df %>% group_by(bin) %>%
    summarise_all(list(sum)) %>%
    mutate(actual = actual/weight,
           pred = pred/weight,
           AvE = actual/pred ) %>%
    arrange(bin)
  return(lift_df_agg)
}

KT_resample_lift = function(n,pred,actual,weight,nbin){
  df = data.frame(actual = actual, weight = weight, predicted =pred)
  lift_df =KT_calc_lift(pred = df$predicted,actual = df$actual,weight = df$weight , nbin=nbin)
  lift_ave = data.frame(bin = lift_df$bin, weight = lift_df$weight ,main =lift_df$AvE )
  lift_actual = data.frame(bin = lift_df$bin, weight = lift_df$weight ,main =lift_df$actual )
  lift_pred = data.frame(bin = lift_df$bin, weight = lift_df$weight ,main =lift_df$pred )
  count = 0
  for (x in sample(1:n*33, n, replace=FALSE) ){
    set.seed(x)
    count=count+1
    test = df[sample(nrow(df), size =  nrow(df) , replace = T), ]
    lift_calc =KT_calc_lift(pred = test$predicted,actual = test$actual,weight = test$weight , nbin=nbin)
    lift_ave[,glue("{count}")] = lift_calc$AvE
    lift_actual[,glue("{count}")] = lift_calc$actual
    lift_pred[,glue("{count}")] = lift_calc$pred
  }
  lift_ave$metric="ave"
  lift_pred$metric = "pred"
  lift_actual$metric = "actual"
  lift = rbind(lift_ave,lift_pred,lift_actual)
  return( lift)
}


KT_plot_lift = function(n , pred, actual ,weight,nbin , title) {
  lift = KT_resample_lift( n = n , pred = pred, actual = actual,nbin = nbin, weight = weight)
  lift = lift %>% 
         rowwise() %>% mutate(lb  =quantile(c_across(3:n+3) , 0.05),
                         ub  =quantile(c_across(3:n+3) , 0.95)) 
  
  ave_plot = lift %>%filter(metric == "ave") %>%
    ggplot(.,aes(x = bin , y = main , group = metric, fill = metric , col = metric )) + geom_line(size = 1.5) +
    geom_point(size = 2.8)+
    geom_ribbon(aes(ymin = lb, ymax = ub ), alpha = 0.3  ,  color = NA)+
    scale_colour_manual("", values = c("red", "blue")) +
    scale_fill_manual("", values = c("red", "blue"))+
    xlab("equal weight banded predictions from low to high")+
    ylab("AvE")+
    ggtitle(title)
  
  lift_plot = lift %>%filter(metric != "ave") %>%
    ggplot(.,aes(x = bin , y = main , group = metric, fill = metric , col = metric)) + geom_line(size = 1.5) +
    geom_point(size = 2.8)+
    geom_ribbon(aes(ymin = lb, ymax = ub ), alpha = 0.3  ,  color = NA)+
    scale_colour_manual("", values = c("red", "blue")) +
    scale_fill_manual("", values = c("red", "blue"))+
    xlab("equal weight banded predictions from low to high")+
    ylab("Predictions and actual")+
    ggtitle(title)
  
  
  return(list(lift_df=  lift, plot= list(lift_plot = lift_plot,ave_plot=ave_plot) ))

}
######################## Double lift ############################

KT_calc_dl = function( actual, weight, base  , challenger , nbin){
  df = data.frame(actual = actual , weight = weight, base = base , challenger = challenger)
  df  %>% 
    filter(weight>0) %>%
    mutate(model_ratio = base/challenger) %>%
    mutate(base=base*weight * (sum(actual)/sum(base*weight)),
           challenger=challenger*weight * (sum(actual)/sum(challenger*weight)))  %>%
    arrange(model_ratio) %>%
    mutate(bin = KT_create_equal_bin(weight,nbin)) %>%
    select(-model_ratio) %>%
    group_by(bin) %>%
    summarise_all(list(sum)) %>%
    mutate(actual = actual/weight,
           base = base/weight,
           challenger = challenger/weight,
           challenger_rb= challenger/base,
           actual_rb = actual / base ) %>%
    arrange(bin)
}

KT_resample_dl = function(n,actual, base  , challenger , weight, nbin ){
  
  dl_sim = list()
  df = data.frame(actual,base,challenger,weight )
  main_dl = KT_calc_dl(
    actual = df$actual,
    base = df$base,
    challenger = df$challenger,
    weight = df$weight,
    nbin = nbin)
  main_dl$sample = "main" 
  dl_sim[["iter_0"]] <- main_dl
  for(x in  seq(1:n)){
    set.seed(x)
    temp = df %>% sample_frac(.,size = 1, replace = T)
    
    dl_sim[[glue("iter_{x}")]] <- data.frame(KT_calc_dl(
      actual = temp$actual,
      base = temp$base,
      challenger = temp$challenger,
      weight = temp$weight,
      nbin = nbin),
      sample = x )
    
  }
  
  variables =  list()
  for (var in c("actual","base", "challenger", "actual_rb", "challenger_rb")){
    
    variables[[var]] <- data.table::rbindlist(dl_sim) %>% 
      select(bin, !!as.name(var),sample) %>%
      pivot_wider(names_from = sample,values_from = !!as.name(var)) %>%
      rowwise() %>%
      mutate(lb  =quantile(c_across(2:n+1) , 0.05, na.rm= T),
             ub  =quantile(c_across(2:n+1) , 0.95 , na.rm= T)) %>%
      select(bin,main,lb,ub) %>%
      mutate(variable = var) %>% 
      ungroup()
    
  }
  return(list(dl_df =  data.frame( data.table::rbindlist(variables),
                                   weight = main_dl$weight) , 
              main_dl =   main_dl)) 
}


KT_plot_dl = function(n, actual, weight, base  , challenger , nbin  ){
  
  test = KT_resample_dl(actual = actual,
                        weight  = weight,
                        challenger =  challenger,
                        base =  base,
                        nbin = nbin,
                        n=n)
  p1 = test$dl_df %>% filter(! grepl("_rb",variable)) %>% 
    mutate(weight= ifelse(variable=="actual",weight,0)) %>%
    ggplot(.,aes(x=bin, group = variable, colour = variable , fill = variable))+
    geom_point(aes(y=main, group = variable), size =2)+
    geom_line(aes(y=main, group = variable), size = 1.2)+
    geom_ribbon(aes(ymin = lb, ymax = ub , group = variable, col = variable , fill = variable ), alpha = 0.3  ,color = NA)+
    ggtitle("Double lift plot")+
    xlab("equal weight banded model ratio from low to high")+
    theme_gray(base_size = 17)+
    theme(
               axis.title.y=element_blank())
    
  
  p2 = test$dl_df %>% filter( grepl("_rb",variable)) %>% 
    mutate(weight= ifelse(variable=="actual",weight,0)) %>%
    ggplot(.,aes(x=bin, group = variable, colour = variable , fill = variable))+
    geom_point(aes(y=main, group = variable), size =2)+
    geom_line(aes(y=main, group = variable), size = 1.2)+
    geom_ribbon(aes(ymin = lb, ymax = ub , group = variable, col = variable , fill = variable ), alpha = 0.3  ,color = NA)+
  ggtitle("Double lift plot - AvE view")+
    scale_colour_manual("", values = c("red", "blue")) +
    scale_fill_manual("", values = c("red", "blue"))+
    geom_hline(yintercept = 1 ,linetype = 2)+
    xlab("equal weight banded model ratio from low to high")+
    theme_gray(base_size = 17)+
    ylab("AvE")
  
  return(list( dl_plot = p1 , dl_rb_plot = p2, dl_df = test ))
}
######################## Triple lift ############################

KT_calc_tl = function( actual, weight, base  , challenger ,challenger2, nbin){
  df = data.frame(actual = actual , weight = weight, base = base , challenger = challenger,challenger2 = challenger2)
  df  %>% 
    filter(weight>0) %>%
    mutate(model_ratio = challenger2/challenger) %>%
    mutate(base=base*weight * (sum(actual)/sum(base*weight)),
           challenger=challenger*weight * (sum(actual)/sum(challenger*weight)),
           challenger2=challenger2*weight * (sum(actual)/sum(challenger2*weight)))  %>%
    arrange(model_ratio) %>%
    mutate(bin = KT_create_equal_bin(weight,nbin)) %>%
    select(-model_ratio) %>%
    group_by(bin) %>%
    summarise_all(list(sum)) %>%
    mutate(actual = actual/weight,
           base = base/weight,
           challenger = challenger/weight,
           challenger_rb= challenger/actual,
           challenger2 = challenger2/weight,
           challenger2_rb= challenger2/actual,
           base_rb = base / actual ) %>%
    arrange(bin)
}

KT_resample_tl = function(n,actual, base  , challenger ,challenger2, weight, nbin ){
  
  tl_sim = list()
  df = data.frame(actual,base,challenger,challenger2,weight )
  main_tl = KT_calc_tl(
    actual = df$actual,
    base = df$base,
    challenger = df$challenger,
    challenger2 = df$challenger2,
    weight = df$weight,
    nbin = nbin)
  main_tl$sample = "main" 
  tl_sim[["iter_0"]] <- main_tl
  for(x in  seq(1:n)){
    set.seed(x)
    temp = df %>% sample_frac(.,size = 1, replace = T)
    
    tl_sim[[glue("iter_{x}")]] <- data.frame(KT_calc_tl(
      actual = temp$actual,
      base = temp$base,
      challenger = temp$challenger,
      challenger2 = temp$challenger2,
      weight = temp$weight,
      nbin = nbin),
      sample = x )
    
  }
  
  variables =  list()
  for (var in c("actual","base", "challenger","challenger2",  "base_rb", "challenger_rb"  , "challenger2_rb" )){
    
    variables[[var]] <- data.table::rbindlist(tl_sim) %>% 
      select(bin, !!as.name(var),sample) %>%
      pivot_wider(names_from = sample,values_from = !!as.name(var)) %>%
      rowwise() %>%
      mutate(lb  =quantile(c_across(2:n+1) , 0.05, na.rm= T),
             ub  =quantile(c_across(2:n+1) , 0.95 , na.rm= T)) %>%
      select(bin,main,lb,ub) %>%
      mutate(variable = var) %>% 
      ungroup()
    
  }
  return(list(tl_df =  data.frame( data.table::rbindlist(variables),
                                   weight = main_tl$weight) , 
              main_tl =   main_tl)) 
}


KT_plot_tl = function(n, actual, weight, base  , challenger ,challenger2  , nbin  ){
  
  test = KT_resample_tl(actual = actual,
                        weight  = weight,
                        challenger =  challenger,
                        challenger2 =  challenger2,
                        base =  base,
                        nbin = nbin,
                        n=n)
  p1 = test$tl_df %>% filter(! grepl("_rb",variable)) %>% 
    mutate(weight= ifelse(variable=="actual",weight,0)) %>%
    ggplot(.,aes(x=bin, group = variable, colour = variable , fill = variable))+
    geom_point(aes(y=main, group = variable), size =2)+
    geom_line(aes(y=main, group = variable), size = 1.2)+
    geom_ribbon(aes(ymin = lb, ymax = ub , group = variable, col = variable , fill = variable ), alpha = 0.2  ,color = NA)+
    ggtitle("Triple lift plot")+
    xlab("equal weight banded model ratio from low to high")+
    theme_gray(base_size = 17)+
    theme(
      axis.title.y=element_blank())
  
  
  p2 = test$tl_df %>% filter( grepl("_rb",variable)) %>% 
    mutate(weight= ifelse(variable=="actual",weight,0)) %>%
    ggplot(.,aes(x=bin, group = variable, colour = variable , fill = variable))+
    geom_point(aes(y=main, group = variable), size =2)+
    geom_line(aes(y=main, group = variable), size = 1.2)+
    geom_ribbon(aes(ymin = lb, ymax = ub , group = variable, col = variable , fill = variable ), alpha = 0.2 ,color = NA)+
    ggtitle("Triple lift plot - AvE view")+
    scale_colour_manual("", values = c("red", "blue" , "green")) +
    scale_fill_manual("", values = c("red", "blue","green"))+
    geom_hline(yintercept = 1 ,linetype = 2)+
    xlab("equal weight banded model ratio from low to high")+
    theme_gray(base_size = 17)+
    ylab("AvE")
  
  return(list( tl_plot = p1 , tl_rb_plot = p2, tl_df = test ))
}

################################ AvE ###################################

KT_calc_ave = function(ft,actual,pred,challenger, weight){
  
  if (missing(challenger)){
    challenger = pred
  }
  pred =   pred*(sum(actual)/sum(pred*weight )) # rebase
  challenger =   challenger*(sum(actual)/sum(challenger*weight )) # rebase
  df = data.frame(ft,actual,pred,challenger,weight )
  df %>% 
    mutate_at(vars(c("pred", "challenger")) , ~.x*weight) %>%
    group_by(ft) %>%
    summarise_all(list(sum)) %>%
    mutate(actual=actual/weight,
           pred=pred/weight,
           challenger=challenger/weight,
           ave = actual/pred,
           challenger_ave = actual/challenger
    ) -> df
  return(df)
  
}
KT_resample_ave = function( n, ft,actual,pred, challenger,weight) {
  if (missing(challenger)){
    challenger = pred
  }
  ave_sim = list()
  df = data.frame(ft,actual,pred,challenger,weight)
  main_ave = KT_calc_ave(ft = df$ft,
                         actual = df$actual,
                         pred = df$pred,
                         challenger = df$challenger,
                         weight = df$weight)
  main_ave$sample = "main" 
  ave_sim[["iter_0"]] <- main_ave
  for(x in  seq(1:n)){
    set.seed(x)
    temp = df %>% sample_frac(.,size = 0.3, replace = FALSE)
    
    ave_sim[[glue("iter_{x}")]] <- data.frame(KT_calc_ave(temp$ft,
                                                          temp$actual,
                                                          temp$pred,
                                                          temp$challenger, 
                                                          temp$weight),
                                              sample = x )
    
  }
  
  variables =  list()
  for (var in c("actual","pred","ave", "challenger", "challenger_ave")){
    
    variables[[var]] <- data.table::rbindlist(ave_sim) %>% 
      select(ft,!!as.name(var),sample) %>%
      pivot_wider(names_from = sample,values_from = !!as.name(var)) %>%
      rowwise() %>%
      mutate(lb  =quantile(c_across(2:n+1) , 0.05, na.rm= T),
             ub  =quantile(c_across(2:n+1) , 0.95 , na.rm= T)) %>%
      select(ft,main,lb,ub) %>%
      mutate(variable = var) %>% 
      ungroup()
    
  }
  
  return( list(ave_df =  data.frame( data.table::rbindlist(variables),
                                     weight = main_ave$weight) , 
               main_ave =   main_ave))
}




KT_plot_ave = function(n, ft,actual,pred, challenger,weight,factor_name,title,rescale=30){
  
  if (missing(challenger)){
    challenger = pred
  }
  
  test = KT_resample_ave( n = n ,
                          ft = ft,
                          actual=actual ,
                          pred = pred,
                          challenger=challenger,
                          weight=weight )
  line_size = 1.2
  point_size = 2.3
  p1 = test$ave_df %>% filter(grepl("actual|pred" , variable)) %>%
    mutate(weight= ifelse(variable=="actual",weight,0)) %>%
    ggplot(.,aes(x=ft, group = variable, colour = variable , fill = variable))+
    geom_point(aes(y=main, group = variable), size =point_size)+
    geom_line(aes(y=main, group = variable), size = line_size)+
    geom_ribbon(aes(ymin = lb, ymax = ub , group = variable, col = variable , fill = variable ), alpha = 0.15  ,color = NA)+
    scale_colour_manual("", values = c("red", "blue")) +
    scale_fill_manual("", values = c("red", "blue"))+
    xlab(factor_name)+
    # ylab("")+
    ggtitle(title)+
    theme_bw() +theme(panel.background = element_blank())+
    # theme_gray(base_size = 17)+
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
    theme(legend.position = "bottom")+
    geom_bar( aes(y=weight/rescale), stat="identity", size=.1, color="black", alpha=.4) +
    scale_y_continuous(name = "",sec.axis = sec_axis(~.*rescale, name="weight"))
  
  p2 = test$ave_df %>% filter(variable %in% c("actual","challenger","pred")) %>%
    mutate(weight= ifelse(variable=="actual",weight,0)) %>%
    ggplot(.,aes(x=ft, group = variable, colour = variable , fill = variable))+
    geom_point(aes(y=main, group = variable), size =point_size)+
    geom_line(aes(y=main, group = variable), size = line_size)+
    geom_ribbon(aes(ymin = lb, ymax = ub , group = variable, col = variable , fill = variable ), alpha = 0.15  ,color = NA)+
    scale_colour_manual("", values = c("red", "green","blue" )) +
    scale_fill_manual("", values = c("red", "green",  "blue" ))+
    xlab(factor_name)+
    # ylab("actual and expected")+
    ggtitle(title)+
    theme_bw() +theme(panel.background = element_blank())+
    # theme_gray(base_size = 17)+
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
    theme(legend.position = "bottom")+
    geom_bar( aes(y=weight/rescale), stat="identity", size=.1, color="black", alpha=.4) +
    scale_y_continuous(name = "",sec.axis = sec_axis(~.*rescale, name="weight"))
  rescale2= rescale*1000
  p3 =test$ave_df %>% filter(grepl("ave" , variable)) %>%
    mutate(weight= ifelse(variable=="ave",weight,0)) %>%
    ggplot(.,aes(x=ft, group = variable, colour = variable , fill = variable))+
    geom_point(aes(y=main, group = variable), size =point_size)+
    geom_line(aes(y=main, group = variable), size = line_size)+
    geom_ribbon(aes(ymin = lb, ymax = ub , group = variable, col = variable , fill = variable ), alpha = 0.15  ,color = NA)+
    scale_colour_manual("", values = c("red", "blue" )) +
    scale_fill_manual("", values = c("red", "blue" ))+
    geom_hline(yintercept = 1 ,linetype = 2)+
    xlab(factor_name)+
    # ylab("actual/expected")+
    ggtitle(title)+
    theme_bw() +theme(panel.background = element_blank())+
    # theme_gray(base_size = 17)+
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
    theme(legend.position = "bottom")+
    geom_bar( aes(y=weight/rescale2), stat="identity", size=.1, color="black", alpha=.4) +
    scale_y_continuous(name = "actual/expected",sec.axis = sec_axis(~.*rescale2, name="weight"))
  
  p4 =test$ave_df %>% filter(variable == "ave") %>%
    mutate(weight= ifelse(variable=="ave",weight,0)) %>%
    ggplot(.,aes(x=ft, group = variable, colour = variable , fill = variable))+
    geom_point(aes(y=main, group = variable), size =point_size)+
    geom_line(aes(y=main, group = variable), size = line_size)+
    geom_ribbon(aes(ymin = lb, ymax = ub , group = variable, col = variable , fill = variable ), alpha = 0.15  ,color = NA)+
    scale_colour_manual("", values = c("red", "blue" )) +
    scale_fill_manual("", values = c("red", "blue" ))+
    geom_hline(yintercept = 1 ,linetype = 2)+
    xlab(factor_name)+
    # ylab("actual/expected")+
    ggtitle(title)+
    theme_bw() +theme(panel.background = element_blank())+
    # theme_gray(base_size = 17)+
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
    theme(legend.position = "bottom")+
    geom_bar( aes(y=weight/rescale2), stat="identity", size=.1, color="black", alpha=.4) +
    scale_y_continuous(name = "actual/expected",sec.axis = sec_axis(~.*rescale2, name="weight"))
  p5 = test$weight %>% ggplot(.,aes(x=ft,y=weight)) + geom_bar(stat= "identity") +
    xlab(factor_name)+
    ggtitle(title)+
    # theme_gray(base_size = 17)+
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))
  
  
  return(list(compare_plot = list(compare_ave_plot = p2,
                                  compare_ave_plot_rb=p3) ,
              model_plot = list(ave_plot = p1,
                            ave_plot_rb = p4),
              
              weight=p5,
              ave_df = test))
  
}

################################# Explain model ####################################
# 
# use_python("C:\\ProgramData\\anaconda3")
# reticulate::py_run_file("H:\\Restricted Share\\DA P&U\\Tech Modelling\\Users\\Khoa\\RPMtools.py") # Compute SHAP and interactions

# Must compute SHAP using the .py file above before plotting them

KT_plot_shap = function(sv, ft,ft_name,excl){
  
  if (! missing(excl)){
    df = data.frame(sv,ft) %>%
      filter(! ft %in% excl)
  }
  else {
    df = data.frame(sv,ft)
  }
  
  if (!is.numeric(ft)){
    df  %>%
      arrange(ft) %>%
      ggplot(.,aes(x=ft, y=sv  )) +geom_point(alpha= 0.3, size = 1.5, colour = "blue"  , fill = "blue")+
      # theme_gray(base_size = 13)+
      theme_bw() +theme(panel.background = element_blank())+
      theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
      xlab(ft_name)+
      ylab("shap_values")+
      ggtitle(glue("{ft_name} SHAP trend")) -> p 
  }else{
    df  %>%
      arrange(ft) %>%
      mutate(lowess = lowess(x = ft,y=sv, f = 1/10)$y)  %>%
      ggplot(.,aes(x=ft, y=sv  )) +geom_point(alpha= 0.3, size = 1.5, colour = "blue"  , fill = "blue")+
      # theme_gray(base_size = 13)+
      theme_bw() +theme(panel.background = element_blank())+
      theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
      geom_line(aes(x=  ft , y = lowess ),lwd = 1)+
      xlab(ft_name)+
      ylab("shap_values")+
      ggtitle(glue("{ft_name} SHAP trend")) -> p 
  }
  
  return(p)
}


KT_plot_shap_w_interaction =  function(sv, ft,ft_name,excl,interaction){
  if (! missing(excl)){
    df = data.frame(sv,ft) %>%
      filter(! ft %in% excl)
  }
  else {
    df = data.frame(sv,ft,interaction)
  }
  df  %>% 
    group_by(interaction)%>%
    arrange(ft) %>%
    ggplot(.,aes(x=ft, y=sv  , colour = interaction )) +geom_point( alpha= 0.3, size = 1.5 )+
    # scale_color_viridis_c()  +
    theme_bw() +theme(panel.background = element_blank())+
    # theme_gray(base_size = 13)+
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
    xlab(ft_name)+
    ylab("shap_values")   ->p
  
  if(is.numeric(interaction))
    p = p+scale_color_viridis_c() 
  return(p)
  
}

KT_plot_compare_shap = function(sv_base,sv_challenger , base_ft, challenger_ft,ft_name){
  
  df_base = data.frame(sv=sv_base,ft = base_ft , scenario = "base")
  df_challenger = data.frame(sv=sv_challenger,ft = challenger_ft,scenario = "challenger")
  df = rbind(df_base,df_challenger)
  df  %>%
    group_by(scenario) %>%
    arrange(ft) %>%
    mutate(lowess = lowess(x = ft,y=sv)$y)  %>%
    ggplot(.,aes(x=ft, y=sv , group = scenario, colour = scenario) ) +
    geom_point(alpha= 0.1, size = 1.5 ,stroke=0.5,shape = 21)+
    # theme_gray(base_size = 13)+
    theme_bw() +theme(panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=0.9))+
    geom_line(aes(x=  ft , y = lowess , group  = scenario , color  = scenario ),lwd = 1,linetype = 1)+
    scale_colour_manual(values=c('blue','red'))+
    xlab(ft_name)+
    ylab("shap_values")+
    ggtitle(glue("{ft_name} SHAP trend"))
}

########################################## UK Map #################################################

KT_prepare_uk_lookup_map = function(){
  postcode_lookup_geometry = list()
  postcode_lookup_poly = list()
  postcode_regex = list()
  postcode_lookup_shp = list()
  lvl = c("area" = "([A-Z][A-Z]{0,1})" ,
          "district"="(([A-Z][A-Z]{0,1})[0-9][A-Z0-9]{0,1})" ,
          "sector"="((([A-Z][A-Z]{0,1})[0-9][A-Z0-9]{0,1}) {0,}[0-9])",
          "postcode"= "^(((([A-Z][A-Z]{0,1})[0-9][A-Z0-9]{0,1}) {0,}[0-9])[A-Z]{2})$")
  postcode = data.table::fread("H:/Restricted Share/DA P&U/Tech Modelling/Users/Khoa/ukpostcodes/ukpostcodes.csv")# we have NI in this file.
  for (dummy in c("area","district","sector" )){
    lookup_poly  = postcode %>% filter(!is.na(longitude)) %>%
      mutate(name =stringr::str_extract(postcode ,lvl[dummy])  ) %>%
      select(name , longitude , latitude) %>%
      group_by(name) %>%
      summarise_all(list(mean)) 
    lookup_geo = lookup_poly %>%
      filter(latitude <90) %>%
      st_as_sf(.,coords = c(2:3))
    postcode_lookup_poly[[dummy]] = lookup_poly
    postcode_lookup_geometry[[dummy]] = lookup_geo
    postcode_regex[[dummy]] = lvl[dummy]
    postcode_lookup_shp[[dummy]] = read_sf(glue("H:\\Restricted Share\\DA P&U\\Tech Modelling\\Users\\Khoa\\UK-postcode-boundaries-Jan-2015\\Distribution\\{dummy}s.shp") )}
    
  return(list(postcode_lookup_poly=postcode_lookup_poly,
              postcode_lookup_geometry=postcode_lookup_geometry,
              postcode_lookup_shp = postcode_lookup_shp,
              postcode_regex=postcode_regex))
}

KT_plot_uk_map = function(df,value, title, size,group = NA, alpha =0.5 ,nrow = 1  ){
  
     point_plot = df %>% ggplot(.) + 
       geom_point(aes(color = !!as.name(value),  
                      geometry = geometry , 
                      size = size) , 
                  alpha = alpha , 
                  stat = "sf_coordinates", 
                  size = size, stroke = 0.15  , 
                  show.legend = T )+
      scale_color_viridis_c(alpha = 1, begin = 0, end = 1)+
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA ),
            legend.position = "bottom",
            plot.background = element_rect(fill = 'black'),
            panel.grid.minor=element_line(colour="black"),
            panel.grid.major=element_line(colour="black"),
            legend.background = element_rect(fill="black"),
            legend.text = element_text(color ="white"),
            legend.title = element_text(color ="white"),
            text = element_text(colour = "white"),
            axis.text = element_text(colour = "white"),
            strip.text = element_text(size = rel(1))) + ylab("Latitude") + xlab("Longtitude") +
       # facet_wrap(.~group,nrow=nrow,strip.position = "top")+
       ggtitle(title)
  
    shape_plot = df %>% ggplot(.) +
      geom_sf(aes(fill = !!as.name(value),  geometry = geometry ) , colour = NA) +
      scale_fill_viridis_c(alpha = 1, begin = 0, end = 1)+
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA ),
            legend.position = "bottom",
            plot.background = element_rect(fill = 'black'),
            panel.grid.minor=element_line(colour="black"),
            panel.grid.major=element_line(colour="black"),
            legend.background = element_rect(fill="black"),
            legend.text = element_text(color ="white"),
            legend.title = element_text(color ="white"),
            text = element_text(colour = "white"),
            axis.text = element_text(colour = "white"),
            strip.text = element_text(size = rel(1))) + ylab("Latitude") + xlab("Longtitude") +
      # facet_wrap(.~group,nrow=nrow,strip.position = "top")+
      ggtitle(title)
    return(list(point_plot = point_plot, shape_plot=shape_plot))
}

################################## XGB Modelling Pipeline #####################################################

# Note working, cores not initialised

# library(doParallel)
# cl <- makeCluster(detectCores())
# registerDoParallel(cl, cores = detectCores())
# clusterExport(cl,c('train_mat' , "test_mat"))
# clusterEvalQ(cl,expr= {
#   library(xgboost)
# })
# 
# stopCluster(cl)
# registerDoSEQ()


KT_xgb_train <- function(train ,
                         train_y ,
                         train_weight ,
                         validate ,
                         validate_y ,
                         validate_weight,
                         params ,
                         verbose = 1) {

  
  train_mat <- xgb.DMatrix(data = as.matrix(train),
                           label = train_y,
                           weight = train_weight)
  
  
  if(missing(validate)){
    watchlist=list(train = train_mat)
    test_mat = NULL
    early_stopping_rounds=NULL
    
  }
  else{
    
    test_mat <- xgb.DMatrix(data = as.matrix(validate),
                            label = validate_y,
                            weight = validate_weight)
    watchlist = list(train = train_mat, test = test_mat)
    early_stopping_rounds = 10
  }
  model <- xgb.train(params = params, 
                     data = train_mat, 
                     nrounds = params$nrounds, 
                     watchlist =watchlist, 
                     print_every_n = 5, 
                     early_stopping_rounds = early_stopping_rounds,
                     maximize = FALSE,
                     verbose = verbose
                     # nthread =5
                     )
  
  if(missing(validate)){
    return(list(model= model))
  }
  else{
    names(model$evaluation_log) <- c("iter" , "train_loss" , "test_loss")
    e <- data.frame(model$evaluation_log)
    e %>% melt(id.vars = "iter") %>% 
      ggplot(.,aes(x=iter,y=value,group = variable ,color = variable)) + geom_line(lwd = 3, alpha = 0.5) -> loss_plot
    # importance <- xgb.importance(model$ft_names, model =model   )
    return(list( model =model ,loss_plot =loss_plot))
  }
}




KT_create_fold_idx <- function(df,k){
  folds = list()
  for(i in 1:k){
    folds[[i]] <-as.integer(seq(i,nrow(df),by = k))
  }
  return(folds)
}

KT_xgb_cv <- function(train,
                      train_y,
                      train_weight,
                      folds ,
                      params,
                      verbose =0){
  
  train_mat <- xgb.DMatrix(data = as.matrix(train),
                           label = train_y,
                           weight = train_weight)
  
  model <- xgb.cv(params = params, 
                  data = train_mat, 
                  nrounds = params$nrounds, 
                  # watchlist = list(train = train_mat, test = test_mat), 
                  print_every_n = 5, 
                  early_stopping_rounds = 10,
                  maximize = FALSE,
                  verbose = verbose,
                  folds = folds)
                  # nfold = length(folds))
  names(model$evaluation_log) <- c("iter" , "train_loss" ,"train_std",  "test_loss", "test_std")
  e <- data.frame(model$evaluation_log)
  e %>% select(iter , train_loss, test_loss) %>% melt(id.vars = "iter") %>% 
  ggplot(.,aes(x=iter,y=value,group = variable ,color = variable)) + geom_line(lwd = 3, alpha = 0.5) -> loss_plot
  return(list( model =model ,loss_plot =loss_plot))
}



KT_xgb_baysian_tune = function(train ,
                               train_y ,
                               train_weight,
                               validate ,
                               validate_y,
                               validate_weight,
                               folds,
                               bounds,
                               nrounds= 400,
                               objective = "reg:tweedie",
                               eval_metric = "tweedie-nloglik@1.5",
                               parallel = F,
                               iters.k = 1,
                               iters.n = 4,
                               initPoints=10,
                               verbose=1){
  
  if(missing(folds)){
    cv = FALSE
  }
  else{
    cv= T
  }
  
  
  
  obj_fun <- function(eta,
                      # nrounds,
                      max_depth, 
                      min_child_weight, 
                      subsample, 
                      colsample_bytree,
                      lambda ,
                      alpha
                      ) {
    
    params = list(eta=eta,
                  nrounds=nrounds,
                  max_depth=max_depth, 
                  min_child_weight=min_child_weight, 
                  subsample=subsample, 
                  colsample_bytree=colsample_bytree, 
                  lambda=lambda,
                  alpha=alpha,
                  # lambda=1,
                  # alpha=0,
                  objective=objective,
                  eval_metric=eval_metric)
    
    if(cv){
      model = KT_xgb_cv(train =train,
                           train_y = train_y,
                           train_weight =train_weight,
                           folds=folds,
                           params = params)$model
      
    }
    else{
      model = KT_xgb_train(train =train,
                           train_y = train_y,
                           train_weight =train_weight,
                           validate  =validate,
                           validate_y  = validate_y,
                           validate_weight = validate_weight,
                           params = params)$model
    }
    
    
    best_iteration = min(which(model$evaluation_log$test_loss == min(model$evaluation_log$test_loss)))
    validate_loss <-    model$evaluation_log[best_iteration, "test_loss"][[1]]
    
    validate_iter = paste( model$evaluation_log[[ "test_loss"]],collapse = "," )
    return(list(Score = as.numeric( validate_loss), num_rounds = best_iteration , validate_iter = validate_iter ))
  }
  
  # Run the Bayesian optimisation
  

  opt_results = bayesOpt(obj_fun, 
                 bounds =bounds , 
                 initPoints = initPoints, 
                 iters.n = iters.n,
                 iters.k = iters.k,
                 parallel =  parallel,
                 verbose = verbose)

  
  tune_iteration = data.frame()
  for (x in  1:nrow(opt_results$scoreSummary)){
    validate_loss = as.numeric(strsplit(opt_results$scoreSummary$validate_iter[x], ",")[[1]])
    train_iteration = seq(1,length(validate_loss))
    BayOpt_iteration=  as.factor(rep(opt_results$scoreSummary[x,"Iteration"][[1]] , length(validate_loss)))
    
    test = data.frame(BayOpt_iteration,train_iteration,validate_loss)
    
    tune_iteration = rbind(tune_iteration,test)
    
  }
  tune_iteration %>% 
    ggplot(.,aes(x=train_iteration , y = validate_loss , colour = BayOpt_iteration, group = BayOpt_iteration))+
    geom_line(lwd = 1.5)  + theme_gray(base_size = 17) -> tune_iteration
  
  hyperparameters = list()
  hyperparameters[["tune_iteration"]] = tune_iteration
  
  # HP = intersect( c("eta","max_depth", "min_child_weight" , "subsample" , "colsample_bytree", "lambda", "alpha")  , colnames(opt_results$scoreSummary)  )
  HP = intersect(names(bounds)  , colnames(opt_results$scoreSummary)  )
  for (x in rev(HP)){
    hyperparameters[[x]]<- opt_results$scoreSummary %>% 
      mutate( Iteration = as.factor( Iteration)) %>%
      ggplot(.,aes(x= !!as.name(x),  y = Score ,  color = Iteration)) + 
      geom_point(size = 2.5) + 
      theme_gray(base_size = 17)
  }
  
  opt_results$scoreSummary = opt_results$scoreSummary %>%  rename(nrounds = num_rounds) 
  opt_results$scoreSummary%>% arrange(Score)  %>% head(1) %>% as.list  -> best_params 
  best_params$objective = objective
  best_params$eval_metric = eval_metric
  scoreSummary <- opt_results$scoreSummary %>% select(c(HP,"Iteration" , "Score","nrounds") ) %>% arrange(Score)
  
  return( list( opt_results=scoreSummary  ,
                hyperparameters_trends=hyperparameters,
                best_params = best_params))
}

KT_xgb_explain = function(model,  pred_data ,excl ,sample_size = 23000){
  
  pred_data = as.matrix(pred_data)
  set.seed(33)
  pred_data_main_effect = pred_data[sample(nrow(pred_data),min(nrow(pred_data) , sample_size), replace = F),]
  shap_main_effect =   predict(model, newdata =pred_data_main_effect, predcontrib = TRUE)  %>% as.data.frame()
  
  for (x in names(shap_main_effect)){
    shap_main_effect[,x] = KT_quantile_clip(shap_main_effect[,x], min = 0.01,max = 0.99)
  }
    
  set.seed(33)
  pred_data_interaction = pred_data[sample(nrow(pred_data),min(nrow(pred_data) , 4000), replace = F),]
  rm(pred_data)
  shap_interaction = data.frame(  predict(model, newdata =pred_data_interaction, predinteraction  = TRUE)) 
  
  return(list(main_effect = list(pred_data_main_effect=pred_data_main_effect,
                                 shap_main_effect=shap_main_effect) ,
              interaction= list(pred_data_interaction=pred_data_interaction,
                                shap_interaction=shap_interaction)
              ))
  
}




#################### Boruta feature selection #####################################

# library(Boruta)
# 
# xgb.boruta=Boruta(train,
#                   y=train_y[[1]],
#                   maxRuns=12, 
#                   doTrace=2,
#                   holdHistory=TRUE,
#                   getImp=getImpXgboost,
#                   max.depth=model$params$max_depth, 
#                   eta=model$params$eta, 
#                   nthread=4, 
#                   min_child_weight=model$params$min_child_weight,
#                   eval_metric=model$params$eval_metric, 
#                   nrounds=model$params$nrounds, 
#                   objective = model$params$objective,
#                   tree_method="hist",
#                   subsample = model$params$subsample,
#                   colsample_bytree = model$params$colsample_bytree,
#                   alpha = model$params$alpha
#                   
#                   
# )
# 
# 
# boruta_dec=attStats(xgb.boruta)
# 
# #get the names of each feature
# imp_features=row.names(boruta_dec)[which(boruta_dec$decision!="Rejected")]
# #get feature importance history
# boruta.imp.df=as.data.frame(xgb.boruta$ImpHistory)
# #keep only confirmed and tentative features
# boruta.imp.df=boruta.imp.df[,names(boruta.imp.df)%in%imp_features]
# #transform the data to a data frame with two columns: feature and importance value
# boruta.imp.df=melt(boruta.imp.df)
# #create a data frame by adding the decision for each feature as well
# boruta.imp.df=cbind.data.frame(boruta.imp.df, 
#                                decision=boruta_dec$decision[match(boruta.imp.df$variable, 
#                                                                   row.names(boruta_dec))])
# #reorder features data frame by the importance median value
# feature_order=with(boruta.imp.df, reorder(variable, value, median, order = TRUE))
# boruta.imp.df$variable=factor(boruta.imp.df$variable, levels = levels(feature_order))
# 
# boruta.imp.df %>% ggplot(.,aes(y = variable ,  x = value , fill =  decision )) + geom_boxplot() + theme(legend.background = "bottom") +theme_gray(base_size = 25)


#############################################################################################################

############# RDR stuff ##########################
# exe_path = "C:/Program Files (x86)/Radar_4_21/RadarCommandLine.exe"
# component_path = "RadarLive_Phoenix1Home.Endcodingmodellingdata.encoded_factor"
# rdr_path = "H:/Restricted Share/DA P&U/Tech Modelling/01 Home/Phase 2/13. R/Peril Name - Radar Home - Phase2 Modelling WorkFlow  v7.rdr"
KT_rdr_cmd= function(exe_path="C:/Program Files (x86)/Radar_4_21/RadarCommandLine.exe", rdr_path, component_path,use_optimiser = F, stdout=F){
  exe_path = glue('&"{exe_path}"')
  rdr_path = glue('"{rdr_path}"')
  component_path = glue('/target:"{component_path}" ')
  license = glue('/emblemLicence:"Standard"')
  edition = glue('/edition:"Optimiser"')
  if(use_optimiser){
    input   = glue("{exe_path} {rdr_path} {component_path} {license} {edition} ")
  }
  else{
    input = glue("{exe_path} {rdr_path} {component_path} {license} ")
  }
  system2("powershell" , 
          input  = input , 
          args = c( "-executionPolicy", "Bypass"),
          stdout=stdout)
}


###################################Useful tools###############################

KT_target_cat_encoding  = function(cat_df , target, weight ){
  if(missing(weight)){
    weight = 1
  }
  dict= list()
  inverse_dict = list()
  df = data.frame(cat_df,target,weight)
  for(x in names(cat_df)){
    df %>% select(target,weight , x) %>%
      group_by(!!as.name(x)) %>%
      summarise_all(list(sum))  %>% 
      mutate(w_avg =target/ weight) %>%
      arrange(w_avg)    %>% 
      mutate(idx = seq(1,nrow(.))) -> agg_df
    dict[[x]] = setNames(agg_df$idx,agg_df[[x]])
    inverse_dict[[x]] = setNames(agg_df[[x]],agg_df$idx)
  }
  
  return(list(dict=dict,
              inverse_dict = inverse_dict))
}

KT_target_cat_encoding_map = function(x, dict ){
  x  %>% map_dbl(~dict[.x])
}

KT_target_cat_encoding_inverse_map = function(x, dict ){
  x %>% map_chr(~dict[.x])
}


KT_get_df_na = function(df){
  df = df %>%
    select_if(function(x) any(is.na(x))) %>%
    summarise_each(funs(sum(is.na(.))/nrow(df)))  %>%
    reshape2::melt() %>%
    arrange(-value)%>%
    rename(NA_proportion = value) %>%
    mutate(NA_count = NA_proportion * nrow(df))
  return(df)
}

KT_plot_df_na = function(df){
  KT_get_df_na(df) %>% 
    ggplot(.,aes(x=reorder(variable,+NA_proportion) ,y= NA_proportion)) +geom_bar(stat = "identity") + 
    coord_flip() + xlab("Features") +
    ggtitle("Missing Proportion")
}


KT_fix_dtype = function(df, date_cols){
  df = df %>% 
    mutate_all(~if_else(.x %in% c("", "null") , NA , .x))
  if( missing(date_cols)){
    df_date = data.frame(matrix(ncol = 0,nrow = nrow(df)))
  }
  
  else {
    df_date = df %>% select(date_cols) %>%
              mutate_all( as.Date ) %>%
             mutate_all( lubridate::ymd )
    df = df %>% select(-date_cols)
  }
  
  df_num = df %>% mutate_all( as.character) %>%
    mutate_all(~if_else(is.na(.x) , "-999999.9999", .x)) %>%
    mutate_all(as.numeric)
  
  char_cols = df_num %>% summarise(across(everything(), ~ sum(is.na(.)))) %>%
    reshape2::melt() %>%
    filter(value >0) %>%
    select(variable) %>%
    pull
  df_num = df_num %>% select(-char_cols) %>%
    mutate_all(~if_else(.x==-999999.9999 ,NA , .x))
  df_char = df %>% select(char_cols) 
  df = cbind(df_num,df_char,df_date)
  return(df)
}

KT_grid_plot  = function(list, ncol){
  do.call("grid.arrange",c(list, ncol=ncol))
}


KT_compare_dist_banded_factors = function(df1,df2,suffix, common_cols) {
  plots = list()
  
  for (ft in common_cols){
    
    p1 =df1 %>% select(!!as.name(ft)) %>%
      
      mutate_all(as.factor) %>%
      count(!!as.name(ft)) %>%
      
      mutate(prop = n/sum(n)) %>%
      select(-n)
    
    p2 = df2 %>% select(!!as.name(ft)) %>%
      mutate_all(as.factor) %>%
      count(!!as.name(ft)) %>%
      mutate(prop = n/sum(n)) %>%
      select(-n)
    
    test = p1 %>% full_join(p2,by =ft , suffix =suffix  )%>%
      reshape2::melt(id.vars =ft ) %>%  
      arrange(!!as.name(ft))
    lvls <- stringr::str_sort(unique(test[,ft]), numeric = TRUE)
    test[,ft] <- factor(test[,ft], levels = lvls)
    
    plots[[ft]]= test %>%  ggplot(.,aes(x=!!as.name(ft) , group = variable, fill = variable, y = value)) + 
      geom_bar(stat= "identity", position = "identity" , alpha = 0.5 ) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
     
  }
  return(plots)
}
KT_clip = function(x,min,max){
  pmin(pmax(x,min),max)
}

KT_quantile_clip = function(x,min,max){
  KT_clip(x,min = quantile(x,min,na.rm= T) , max =  quantile(x,max,na.rm = T))
}

KT_hex_scatter_plot = function(df,x , y, title){
  breaks = c(0.00001,0.0001,0.001,	0.01,		0.1)
  df %>%
    ggplot(.,aes(x = !!as.name(x) ,y = !!as.name(y) )) + 
    
    stat_binhex(aes(fill= {print(sum(..count..));..count../sum(..count..)}),
                colour = "white",
                stat = "identity") +
    
    scale_fill_gradient(name = "proportion of holdout",
                        trans = "log",
                        breaks = breaks, 
                        labels = breaks,
                        low = "lightgrey",
                        
                        high = "darkred")+
    ggtitle(title)}



KT_get_unique_values_in_each_col = function(df, max_level = 10){
  df = data.frame(df)
  rand_sample = data.frame()
  head_sample=data.frame()
  tail_sample=data.frame()
  for (x in colnames(df)){
    if(is.numeric(df[,x])){
      values = round(  unique( df[,x]),3)
    }
    else{
      values = unique(  df[,x])
    }
    n_levels = length(values)
    
    if (n_levels > max_level){
      values_rand_sample = sample(   values , size =max_level , replace = F)
    }
    else {
      values_rand_sample = values
    }
    sorted_values=sort(values, decreasing = T)
    
    values_head_sample = head(sorted_values,max_level)
    values_tail_sample = tail(sorted_values,max_level)
    
    temp_rand_sample  <- data.frame(variable = x,dtype = class(df[,x]) , n_levels = n_levels,  level_name = paste( sort (values_rand_sample,decreasing = T),collapse = "," ))
    rand_sample <-rbind(rand_sample,temp_rand_sample)
    
    temp_head_sample  <- data.frame(variable = x,dtype = class(df[,x]) , n_levels = n_levels,  level_name = paste( values_head_sample,collapse = "," ))
    head_sample <-rbind(head_sample,temp_head_sample)
      
    temp_tail_sample  <- data.frame(variable = x,dtype = class(df[,x]) , n_levels = n_levels,  level_name = paste( values_tail_sample,collapse = "," ))
    tail_sample <-rbind(tail_sample,temp_tail_sample)

  }
  return(list(rand_sample=rand_sample,head_sample=head_sample, tail_sample=tail_sample))

}


KT_plot_factor_dist = function(df ,cols_to_plot ){
  pb = txtProgressBar(min = 0, max = length(cols_to_plot), initial = 0 , style = 3)
 
  df = data.frame(df)
  dist_list = list()
  i = 0
  for (x in  cols_to_plot ){
    i= i+1
    if (!is.numeric(df[,x])){
      lvls <- stringr::str_sort(unique(df[,x]), numeric = TRUE)
      df[,x] <- factor(df[,x], levels = lvls)
      p = df %>% count(!!as.name(x)) %>%
        mutate(prop= n/sum(n)) %>% 
        ggplot(.,aes(x=!!as.name(x), y = prop))+geom_bar(stat = "identity")+ 
        # theme_gray(base_size = 17)
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
      
    }
    else {
      p =  df %>% select(x) %>%
        ggplot(.,aes(x=!!as.name(x)))+ geom_histogram(aes(y=..count../sum(..count..)))+ 
        # theme_gray(base_size = 17)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
        
    }
    
    
    dist_list[[x]] <-p 
    setTxtProgressBar(pb,i)
  }
  return(dist_list)
}


## tab plot example 
# AVE on unseen data  {.tabset .tabset-pills}
# ```{r , include=F}
# plotly::plot_ly()
# ```
# 
# ```{r Plot AvE, results="asis", fig.width=15,fig.height=12}
# 
# 
# 
# for (i in names(AvE_plots)) {
#   cat("## ",i,"\n")
#   # print(htmltools::tagList(plotly::as.widget(plotly::subplot( plotly::ggplotly( AvE_plots[[ft]]$compare_ave_plot )))))
#   print(AvE_plots[[i]]$compare_ave_plot_rb)
#   cat("\n\n")
# }
# ```

# pb = txtProgressBar(min = 0, max = n, initial = 0)
#setTxtProgressBar(pb,iter_i)
# rename_with(~paste0(.x , "ClaimInd"), matches("c$|b$"))

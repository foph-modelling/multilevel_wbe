#:::::::::::::::::::::::::::::
# Project: multilevel_wbe
# description: create figure for deviation from average trend by ARA
# creation: jriou 
# init date: 2023-06-14
#:::::::::::::::::::::::::::::

mw_131_map_deviation_from_average = function(model,corr,ww,shp,ntop=NULL) {
  
  model = ma5.3.3
  corr = corr_all_ara
  shp = shapes
  ww = ww_all
  
  # extract ARA trends
  nara = length(unique(ww$ara1))
  corr_periods = ww %>% 
    dplyr::select(date,period) %>% 
    dplyr::distinct()
  corr_days = tibble(day=0:max(ww$day),
         date=seq.Date(from=min(ww$date),to=max(ww$date),by=1)) %>% 
    left_join(corr_periods,by = join_by(date))
  tt = model$summary.random$day1
  alldays = unique(tt$ID)
  ndays = length(alldays)
  
  tt = tt %>%  
      dplyr::bind_cols(day=rep(alldays,nara),
                     ara1=rep(1:nara,each=ndays)) %>% 
    dplyr::left_join(corr,by = join_by(ara1)) %>% 
    dplyr::left_join(corr_days,by="day") %>% 
    as_tibble()
  tt2 = tt
  # filter
  if(!is.null(ntop)) {
    extremes = tt %>% 
      dplyr::filter(`0.025quant`>0 | `0.975quant`<0) %>% 
      dplyr::group_by(ara_name) %>% 
      dplyr::summarise(max=max(mean),
                       min=min(mean),
                       abs=max(max,-min),
                       .groups="keep") %>%
      ungroup() %>% 
      dplyr::mutate(overrank=min_rank(-abs)) %>% 
      dplyr::filter(overrank<=ntop)  %>% 
      arrange(-abs) 
     
    tt2 = dplyr::filter(tt,ara_name %in% extremes$ara_name) %>% 
      dplyr::mutate(ara_name=factor(ara_name,levels=extremes$ara_name)) 
  }
  g1 = tt2 %>% 
    ggplot() +
    geom_hline(yintercept=1,linetype=2,alpha=.5) +
    geom_ribbon(aes(x=date,ymin=exp(`0.025quant`),ymax=exp(`0.975quant`),fill=NUTS2_name),alpha=.5) +
    geom_line(aes(x=date,y=exp(mean),colour=NUTS2_name)) +
    facet_wrap(~ara_name) +
    scale_colour_discrete(guide="none") +
    scale_fill_discrete(guide="none") +
    scale_y_continuous(trans="log",breaks = c(.1,1,10)) +
    coord_cartesian(ylim=c(.05,20)) +
    labs(title=paste0("Top deviations from average time trend (top ",ntop,")"),x="Day",y="Relative viral load by ARA") + 
    theme(axis.text.x = element_text(angle=45,hjust=1))
  # print(g1)
  # map max
  tt3 = tt %>% 
    dplyr::filter(`0.025quant`>0) %>%
    dplyr::group_by(ara_name,ara_id) %>% 
    dplyr::summarise(value=max(mean),
                     .groups="keep") %>% 
    ungroup()
  mm = shp$ara_shp %>% 
    dplyr::left_join(tt3, by = join_by(ara_id)) %>% 
    mutate(var="Above")
  # map min
  tt4 = tt %>% 
    dplyr::filter(`0.975quant`<0) %>%
    dplyr::group_by(ara_name,ara_id) %>% 
    dplyr::summarise(value=min(mean),
                     .groups="keep") %>% 
    ungroup()
  nn = shp$ara_shp %>% 
    dplyr::left_join(tt4, by = join_by(ara_id)) %>% 
    mutate(var="Below")
  # plot max min deviation
  pp = bind_rows(mm,nn)
  g2 = ggplot() +
    geom_sf(data=shp$canton_shp,fill="grey95",colour="grey70") +
    geom_sf(data=shp$see_shp,fill="white") +
    geom_sf(data=pp,colour="black",aes(fill=exp(value))) +
    facet_wrap(~var,ncol=2) +
    scale_fill_viridis_c(trans="log10",limits=c(0.1,10)) +
    labs(title="Maximal and minimal deviations from average time trend",fill="Relative VL") 
  # print(g2)
  
  # map period
  qq = tt %>% 
    dplyr::filter(!is.na(period)) %>% 
    dplyr::group_by(ara_name,ara_id,period) %>% 
    dplyr::summarise(mean=mean(mean),
                     start=min(date),
                     end=max(date),
                     .groups="drop",)
  oo = shp$ara_shp %>% 
    dplyr::left_join(qq, by = join_by(ara_id)) %>% 
    mutate(var="Below")
  g3 = ggplot() +
    geom_sf(data=shp$canton_shp,fill="grey95",colour="grey70") +
    geom_sf(data=shp$see_shp,fill="white") +
    geom_sf(data=oo,colour="black",aes(fill=exp(mean))) +
    facet_wrap(~period) +
    scale_fill_viridis_c(trans="log10",limits=c(0.1,10)) +
    labs(title="Maximal and minimal deviations from average time trend",fill="Relative VL") 
  return(list(g1,g2,g3))

}
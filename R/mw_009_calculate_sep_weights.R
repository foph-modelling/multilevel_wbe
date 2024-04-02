
if(FALSE) {
swissboundaries_BFS_NATION = st_read('data/spatial/bfs/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_LANDESGEBIET.shp')
bfs_crs = st_crs(swissboundaries_BFS_NATION)
pop_data = fread('data/population_statistics/STATPOP2021.csv')
pop_data_slim = pop_data[,c('E_KOORD', 'N_KOORD',  'B21BTOT')]

pop_sf = st_as_sf(pop_data_slim, coords = c("E_KOORD", "N_KOORD"), 
                            crs = bfs_crs, remove=FALSE) 


sep_sf = readRDS("../../02_data/sep3/ssep3_user_geo.Rds") %>% 
  sf::st_transform(crs=bfs_crs)


joint_sf <- st_join(sep_sf, pop_sf, join = st_nearest_feature)

joint_dt = data.table(st_drop_geometry(joint_sf))
head(joint_dt)
joint_dt[, hect_id := paste0(E_KOORD, N_KOORD)]
joint_dt[, sep_pop_weight := mean(B21BTOT)/.N, by=c('hect_id')]

write.csv(joint_dt[, c('gisid', 'sep_pop_weight')], file = 'data/sep_weights.csv')
}
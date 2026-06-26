library(ggrepel)

yx_cents <- yx_mcda |>
  st_centroid() |>
  bind_cols(st_coordinates(st_centroid(yx_mcda)) |> 
              as.data.frame() |> 
              rename(lon = X, lat = Y))

dl_cents <- dl_mcda |>
  st_centroid() |>
  bind_cols(st_coordinates(st_centroid(dl_mcda)) |> 
              as.data.frame() |> 
              rename(lon = X, lat = Y))


yuexiu_names <- c(
  "滨江街道" = "Binjiang", "矿泉街道" = "Kuangquan",
  "珠光街道" = "Zhuguang", "建设街道" = "Jianshe",
  "光塔街道" = "Guangta",  "大东街道" = "Dadong",
  "猎德街道" = "Liede",    "农林街道" = "Nonglin",
  "华乐街道" = "Huayue",   "登峰街道" = "Dengfeng",
  "洪桥街道" = "Hongqiao", "六榕街道" = "Liurong",
  "流花街道" = "Liuhua",   "白云街道" = "Baiyun",
  "人民街道" = "Renmin",   "北京街道" = "Beijing Jie",
  "东山街道" = "Dongshan", "梅花村街道" = "Meihuacun",
  "黄花岗街道" = "Huanghuagang", "桂花岗街道" = "Guihuagang",
  "诗书街道" = "Shishu",   "寺右街道" = "Siyou",
  "素社街道" = "Sushe",    "五仙观街道" = "Wuxianguan",
  "西湖路街道" = "Xihu Lu","文德路街道" = "Wende Lu",
  "华林街道" = "Hualin",
  # newly added
  "大塘街道"   = "Datang",
  "龙津街道"   = "Longjin",
  "金花街道"   = "Jinhua",
  "站前街道"   = "Zhanqian",
  "沙河街道"   = "Shahe",
  "沙东街道"   = "Shadong",
  "三元里街道" = "Sanyuanli",
  "景泰街道"   = "Jingtai"
)

yx_cents <- yx_cents |>
  mutate(name_en = yuexiu_names[name],
         name_en = ifelse(is.na(name_en), name, name_en))


p_label_yx <- ggplot(yx_mcda) +
  geom_sf(aes(fill = priority_tier), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = tier_colours, name = "NbS Priority") +
  geom_label_repel(data = yx_cents,
                   aes(x = lon, y = lat, label = name_en),
                   size = 2.5, max.overlaps = 30,
                   box.padding = 0.3, label.size = 0.1,
                   fill = "white", alpha = 0.85) +
  theme_minimal(base_size = 10) +
  labs(title = "Yuexiu Subdistricts — Named Reference Map",
       subtitle = "Coloured by NbS priority tier")

p_label_dl <- ggplot(dl_mcda) +
  geom_sf(aes(fill = priority_tier), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = tier_colours, name = "NbS Priority") +
  geom_label_repel(data = dl_cents,
                   aes(x = lon, y = lat, label = wijknaam),
                   size = 2.5, max.overlaps = 20,
                   box.padding = 0.3, label.size = 0.1,
                   fill = "white", alpha = 0.85) +
  theme_minimal(base_size = 10) +
  labs(title = "Delft Wijken — Named Reference Map",
       subtitle = "Coloured by NbS priority tier")

ggsave(file.path(OUT_ROOT, "fig_reference_maps.png"),
       p_label_yx + p_label_dl, width = 16, height = 8, dpi = 300)
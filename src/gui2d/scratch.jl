using Pkg
pkg"activate ."
using NMRAnalysis
using GLMakie

ENV["JULIA_DEBUG"]=NMRAnalysis

##
relaxation2d("/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11",
             "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/lists/vd/t2-vd.cw")

# t = [0.01568
#     0.03136
#     0.04704
#     0.06272
#     0.0784
#     0.09408
#     0.12544
#     0.1568]
# relaxation2d("/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11", t)

# relaxation2d(
#     ["/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11",
#     "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/12"],
#     ["/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/lists/vd/t2-vd.cw",
#     "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/12/lists/vd/t2-vd.cw"])

# relaxation2d(
#     ["/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/pdata/231",
#      "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/pdata/232",
#      "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/pdata/233",
#      "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/pdata/234",
#      "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/pdata/235",
#      "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/pdata/236",
#      "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/pdata/237",
#      "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/pdata/238",
#     ], "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/11/lists/vd/t2-vd.cw")

# relaxation2d("/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/17",
#     "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/17/lists/vd/t1.cw")

##
hetnoe2d(["/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/21/pdata/231",
          "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/21/pdata/232",
          "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/23/pdata/231",
          "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/23/pdata/232",
          "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/25/pdata/231",
          "/Users/chris/NMR/crick-950/kleo_CRT_CTD_relax_241218/25/pdata/232"],
         [false, true, false, true, false, true])

##
pre2d(["/Users/chris/NMR/crick-800/chris_hewl_231106/18/",
       "/Users/chris/NMR/crick-800/chris_hewl_231106/44/",
       "/Users/chris/NMR/crick-800/chris_hewl_231106/34/",
       "/Users/chris/NMR/crick-800/chris_hewl_231106/24/"],
      [0.0, 1.0, 2.0, 5.0],
      :hmqc,
      0.015)

renderer=bin/path6d
png=bin/hdrtopng
grid=bin/grid
stats=bin/stats

mkdir cornell_readme

## 2d
# random
$renderer --prefix cornell_readme/cornell-random --sampler random -s   1 -d 7 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell-random --sampler random -s   4 -d 7 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell-random --sampler random -s  16 -d 7 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell-random --sampler random -s  64 -d 7 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell-random --sampler random -s 256 -d 7 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help

# sobol++
$renderer --prefix cornell_readme/cornell-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s   1 -d 2 --remap_dims "7 5:0 6:1" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s   4 -d 2 --remap_dims "7 5:0 6:1" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s  16 -d 2 --remap_dims "7 5:0 6:1" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s  64 -d 2 --remap_dims "7 5:0 6:1" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s 256 -d 2 --remap_dims "7 5:0 6:1" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help

## reference 4k
$renderer --prefix cornell_readme/cornell-ref --sampler owen -s 4096 -d 2 --remap_dims "7 5:0 6:1" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --force_direct --scene scenes/cornell.txt --help

## png
$png cornell_readme/cornell-ref-04096.pfm
$png --range cornell_readme/cornell-ref-04096.pfm cornell_readme/cornell-random-0????.pfm
$png --range cornell_readme/cornell-ref-04096.pfm cornell_readme/cornell-sobolpp-0????.pfm


## results
$grid cornell_readme/cornell.png \
	--line 1 --title "ref 4k" cornell_readme/cornell-ref-04096.png \
	--line 2 --title "random" cornell_readme/cornell-random-0????.png \
	--line 3 --title "sobol++" cornell_readme/cornell-sobolpp-0????.png

## errors
$stats cornell_readme/cornell-random.pfm cornell_readme/cornell-ref-04096.pfm cornell_readme/cornell-random-0????.pfm 
$stats cornell_readme/cornell-sobolpp.pfm cornell_readme/cornell-ref-04096.pfm cornell_readme/cornell-sobolpp-0????.pfm 

$grid cornell_readme/errors_cornell.png \
	--line 1 --title "ref 4k" cornell_readme/cornell-ref-04096.png \
	--line 2 --title "random" cornell_readme/cornell-random-0????.png \
	--line 3 cornell_readme/error_cornell-random-0????.png \
	--line 4 --title "sobol++" cornell_readme/cornell-sobolpp-0????.png \
	--line 5 cornell_readme/error_cornell-sobolpp-0????.png

## pixel errors
$stats --pixel 100 100 cornell_readme/cornell-random.pfm cornell_readme/cornell-ref-04096.pfm cornell_readme/cornell-random-0????.pfm 
$stats --pixel 100 100 cornell_readme/cornell-sobolpp.pfm cornell_readme/cornell-ref-04096.pfm cornell_readme/cornell-sobolpp-0????.pfm 

## visualize error plots with gnuplot
gnuplot <<EOF
set term pdfcairo
set output "cornell_readme/pixel_100_100_errors.pdf"
set logscale xy
set xtics 4
set grid
plot "./cornell_readme/pixel_100_100_cornell-random.txt" w lines lw 1.5 t "random", "./cornell_readme/pixel_100_100_cornell-sobolpp.txt" w lines lw 1.5 t "sobol++" 
EOF


## 4d
$renderer --prefix cornell_readme/cornell1-random --sampler random -s   1 -d 9 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell1-random --sampler random -s   4 -d 9 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell1-random --sampler random -s  16 -d 9 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell1-random --sampler random -s  64 -d 9 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell1-random --sampler random -s 256 -d 9 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help

# sobol++
$renderer --prefix cornell_readme/cornell1-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s   1 -d 4 --remap_dims "9 5:0 6:1 7:2 8:3" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell1-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s   4 -d 4 --remap_dims "9 5:0 6:1 7:2 8:3" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell1-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s  16 -d 4 --remap_dims "9 5:0 6:1 7:2 8:3" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell1-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s  64 -d 4 --remap_dims "9 5:0 6:1 7:2 8:3" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell1-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s 256 -d 4 --remap_dims "9 5:0 6:1 7:2 8:3" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help

## reference 4k
$renderer --prefix cornell_readme/cornell1-ref --sampler owen -s 4096 -d 4 --remap_dims "9 5:0 6:1 7:2 8:3" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help

## png
$png cornell_readme/cornell1-ref-04096.pfm
$png --range cornell_readme/cornell1-ref-04096.pfm cornell_readme/cornell1-random-0????.pfm
$png --range cornell_readme/cornell1-ref-04096.pfm cornell_readme/cornell1-sobolpp-0????.pfm


## results
$grid cornell_readme/cornell1.png \
	--line 1 --title "ref 4k" cornell_readme/cornell1-ref-04096.png \
	--line 2 --title "random" cornell_readme/cornell1-random-0????.png \
	--line 3 --title "sobol++" cornell_readme/cornell1-sobolpp-0????.png

## errors
$stats cornell_readme/cornell1-random.pfm cornell_readme/cornell1-ref-04096.pfm cornell_readme/cornell1-random-0????.pfm 
$stats cornell_readme/cornell1-sobolpp.pfm cornell_readme/cornell1-ref-04096.pfm cornell_readme/cornell1-sobolpp-0????.pfm 

$grid cornell_readme/errors_cornell1.png \
	--line 1 --title "ref 4k" cornell_readme/cornell1-ref-04096.png \
	--line 2 --title "random" cornell_readme/cornell1-random-0????.png \
	--line 3 cornell_readme/error_cornell1-random-0????.png \
	--line 4 --title "sobol++" cornell_readme/cornell1-sobolpp-0????.png \
	--line 5 cornell_readme/error_cornell1-sobolpp-0????.png

## pixel errors
$stats --pixel 100 100 cornell_readme/cornell1-random.pfm cornell_readme/cornell1-ref-04096.pfm cornell_readme/cornell1-random-0????.pfm 
$stats --pixel 100 100 cornell_readme/cornell1-sobolpp.pfm cornell_readme/cornell1-ref-04096.pfm cornell_readme/cornell1-sobolpp-0????.pfm 

## visualize error plots with gnuplot
gnuplot <<EOF
set term pdfcairo
set output "cornell_readme/pixel_100_100_errors1.pdf"
set logscale xy
set xtics 4
set grid
plot "./cornell_readme/pixel_100_100_cornell1-random.txt" w lines lw 1.5 t "random", "./cornell_readme/pixel_100_100_cornell1-sobolpp.txt" w lines lw 1.5 t "sobol++" 
EOF


## 6d
$renderer --prefix cornell_readme/cornell2-random --sampler random -s   1 -d 11 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell2-random --sampler random -s   4 -d 11 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell2-random --sampler random -s  16 -d 11 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell2-random --sampler random -s  64 -d 11 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell2-random --sampler random -s 256 -d 11 --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help

# sobol++
$renderer --prefix cornell_readme/cornell2-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s   1 -d 6 --remap_dims "11 5:0 6:1 7:2 8:3 9:4 10:5" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell2-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s   4 -d 6 --remap_dims "11 5:0 6:1 7:2 8:3 9:4 10:5" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell2-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s  16 -d 6 --remap_dims "11 5:0 6:1 7:2 8:3 9:4 10:5" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell2-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s  64 -d 6 --remap_dims "11 5:0 6:1 7:2 8:3 9:4 10:5" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help
$renderer --prefix cornell_readme/cornell2-sobolpp --sampler sobolpp --sobol_table ../../data/sobol_init_tab.dat -s 256 -d 6 --remap_dims "11 5:0 6:1 7:2 8:3 9:4 10:5" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 2 --scene scenes/cornell.txt --help

## reference 4k
$renderer --prefix cornell_readme/cornell2-ref --sampler owen -s 4096 -d 6 --remap_dims "11 5:0 6:1 7:2 8:3 9:4 10:5" --remap_bounce true --source_2d true --brdf_2d true --source_quad1 true --path_length 1 --scene scenes/cornell.txt --help


## png
$png cornell_readme/cornell2-ref-04096.pfm
$png --range cornell_readme/cornell2-ref-04096.pfm cornell_readme/cornell2-random-0????.pfm
$png --range cornell_readme/cornell2-ref-04096.pfm cornell_readme/cornell2-sobolpp-0????.pfm


## results
$grid cornell_readme/cornell2.png \
	--line 1 --title "ref 4k" cornell_readme/cornell2-ref-04096.png \
	--line 2 --title "random" cornell_readme/cornell2-random-0????.png \
	--line 3 --title "sobol++" cornell_readme/cornell2-sobolpp-0????.png

## errors
$stats cornell_readme/cornell2-random.pfm cornell_readme/cornell2-ref-04096.pfm cornell_readme/cornell2-random-0????.pfm 
$stats cornell_readme/cornell2-sobolpp.pfm cornell_readme/cornell2-ref-04096.pfm cornell_readme/cornell2-sobolpp-0????.pfm 

$grid cornell_readme/errors_cornell2.png \
	--line 1 --title "ref 4k" cornell_readme/cornell2-ref-04096.png \
	--line 2 --title "random" cornell_readme/cornell2-random-0????.png \
	--line 3 cornell_readme/error_cornell2-random-0????.png \
	--line 4 --title "sobol++" cornell_readme/cornell2-sobolpp-0????.png \
	--line 5 cornell_readme/error_cornell2-sobolpp-0????.png

## pixel errors
$stats --pixel 100 100 cornell_readme/cornell2-random.pfm cornell_readme/cornell2-ref-04096.pfm cornell_readme/cornell2-random-0????.pfm 
$stats --pixel 100 100 cornell_readme/cornell2-sobolpp.pfm cornell_readme/cornell2-ref-04096.pfm cornell_readme/cornell2-sobolpp-0????.pfm 

## visualize error plots with gnuplot
gnuplot <<EOF
set term pdfcairo
set output "cornell_readme/pixel_100_100_errors2.pdf"
set logscale xy
set xtics 4
set grid
plot "./cornell_readme/pixel_100_100_cornell2-random.txt" w lines lw 1.5 t "random", "./cornell_readme/pixel_100_100_cornell2-sobolpp.txt" w lines lw 1.5 t "sobol++" 
EOF


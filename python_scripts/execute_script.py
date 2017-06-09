from os import system

string = 'python prec_stats.py --date_start 2016052800 --date_end 2016060800 --nens 2 --cld_size_bin_triplet 0 3e9 40 --cld_size_sep_bin_triplet 0 2e9 40 --cld_prec_bin_triplet 0 7e9 40 --cld_prec_sep_bin_triplet 0 4e9 40 --pp_name {} --plot_name {} --cld_y_type relative_frequency --sep_perimeter {} --rdf_cov_thresh 0.01 --time_start 6 --rdf_curve_times 9 15 --no_det --rdf_dr {}'
system(string.format('all-days_2-mems_perim-3_dr-1',
                    'all-days_2-mems_perim-3_dr-1',
                    '3', '1'))
system(string.format('all-days_2-mems_perim-5_dr-1',
                    'all-days_2-mems_perim-5_dr-1',
                    '5', '1'))
system(string.format('all-days_2-mems_perim-11_dr-1',
                    'all-days_2-mems_perim-11_dr-1',
                    '11', '1'))

system(string.format('all-days_2-mems_perim-3_dr-2',
                    'all-days_2-mems_perim-3_dr-2',
                    '3', '2') + ' --which_plot rdf')
system(string.format('all-days_2-mems_perim-5_dr-2',
                    'all-days_2-mems_perim-5_dr-2',
                    '5', '2') + ' --which_plot rdf')
system(string.format('all-days_2-mems_perim-11_dr-2',
                    'all-days_2-mems_perim-11_dr-2',
                    '11', '2') + ' --which_plot rdf')

system(string.format('all-days_2-mems_perim-3_dr-1',
                    'all-days_2-mems_perim-3_dr-1_sep',
                    '3', '1') + ' --which_plot rdf --rdf_sep')
system(string.format('all-days_2-mems_perim-5_dr-1',
                    'all-days_2-mems_perim-5_dr-1_sep',
                    '5', '1') + ' --which_plot rdf')
system(string.format('all-days_2-mems_perim-11_dr-1',
                    'all-days_2-mems_perim-11_dr-1_sep',
                    '11', '1') + ' --which_plot rdf --rdf_sep')

system(string.format('all-days_2-mems_perim-3_dr-2',
                    'all-days_2-mems_perim-3_dr-2_sep',
                    '3', '2') + ' --which_plot rdf --rdf_sep')
system(string.format('all-days_2-mems_perim-5_dr-2',
                    'all-days_2-mems_perim-5_dr-2_sep',
                    '5', '2') + ' --which_plot rdf --rdf_sep')
system(string.format('all-days_2-mems_perim-11_dr-2',
                    'all-days_2-mems_perim-11_dr-2_sep',
                    '11', '2') + ' --which_plot rdf --rdf_sep')
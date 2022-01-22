# prelim_params.py


# set up the landscape
import numpy as np
b4 = np.vstack([np.linspace(1, 0, 50) for _ in range(50)])
af = np.vstack([np.linspace(1, 0.5, 50) for _ in range(50)])
stable = np.vstack([np.linspace(1, 0, 50) for _ in range(50)])
K = np.ones((50,50))

# show the landscape, for debugging, if requested
debug_landscape = False
if debug_landscape:
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    im1 = ax1.imshow(b4, cmap='spring', vmin=0, vmax=1)
    plt.colorbar(im1)
    ax3 = fig.add_subplot(223)
    im3 = ax3.imshow(af, cmap='spring', vmin=0, vmax=1)
    plt.colorbar(im3)
    ax2 = fig.add_subplot(222)
    im2 = ax2.imshow(stable, cmap='winter', vmin=0, vmax=1)
    plt.colorbar(im2)
    ax4 = fig.add_subplot(224)
    im4 = ax4.imshow(K, cmap='autumn', vmin=0, vmax=1)
    plt.colorbar(im4)
    plt.show()


# This is a parameters file generated by Geonomics
# (by the gnx.make_parameters_file() function).


                   #   ::::::          :::    :: :::::::::::#
             #::::::    ::::   :::      ::    :: :: ::::::::::: ::#
          #:::::::::     ::            ::   :::::::::::::::::::::::::#
        #::::::::::                      :::::::::: :::::: ::::::::  ::#
      #  : ::::  ::                    ::::  : ::    :::::::: : ::  :    #
     # GGGGG :EEEE: OOOOO   NN   NN   OOOOO   MM   MM IIIIII  CCCCC SSSSS #
    # GG     EE    OO   OO  NNN  NN  OO   OO  MM   MM   II   CC     SS     #
    # GG     EE   OO     OO NN N NN OO     OO MMM MMM   II   CC     SSSSSS #
    # GG GGG EEEE OO     OO NN  NNN OO     OO MM M MM   II   CC         SS #
    # GG   G EE    OO   OO  NN   NN  OO   OO  MM   MM   II   CC        SSS #
     # GGGGG :EEEE: OOOOO   NN   NN   OOOOO   MM   MM IIIIII  CCCCC SSSSS #
      #    : ::::::::               :::::::::: ::              ::  :   : #
        #:    :::::                    :::::: :::             :::::::  #
          #    :::                      :::::  ::              ::::: #
             #  ::                      ::::                      #
                   #                                        #
                      #  :: ::    :::             #


params = {
###############################################################################

###################
#### LANDSCAPE ####
###################
    'landscape': {

    ##############
    #### main ####
    ##############
        'main': {
            #x,y (a.k.a. j,i) dimensions of the Landscape
            'dim':                      (50,50),
            #x,y resolution of the Landscape
            'res':                      (1,1),
            #x,y coords of upper-left corner of the Landscape
            'ulc':                      (0,0),
            #projection of the Landscape
            'prj':                      None,
            }, # <END> 'main'

    ################
    #### layers ####
    ################
        'layers': {

            #layer name (LAYER NAMES MUST BE UNIQUE!)
            'shift': {

        #-------------------------------------#
        #--- layer num. 0: init parameters ---#
        #-------------------------------------#

                #initiating parameters for this layer
                'init': {

                    #parameters for a 'defined'-type Layer
                    'defined': {
                        #raster to use for the Layer
                        'rast':                   b4,
                        #point coordinates
                        'pts':                    None,
                        #point values
                        'vals':                   None,
                        #interpolation method {None, 'linear', 'cubic',
                        #'nearest'}
                        'interp_method':          None,

                        }, # <END> 'defined'

                    }, # <END> 'init'

            #---------------------------------------#
            #--- layer num. 0: change parameters ---#
            #---------------------------------------#

                #landscape-change events for this Layer
                'change': {

                    0: {
                        #array or file for final raster of event, or directory
                        #of files for each stepwise change in event
                        'change_rast':              af,
                        #starting timestep of event
                        'start_t':          1500,
                        #ending timestep of event
                        'end_t':            1750,
                        #number of stepwise changes in event
                        'n_steps':          250,
                        }, # <END> event 0

                    }, # <END> 'change'

                }, # <END> layer num. 0


            #layer name (LAYER NAMES MUST BE UNIQUE!)
            'stable': {

        #-------------------------------------#
        #--- layer num. 1: init parameters ---#
        #-------------------------------------#

                #initiating parameters for this layer
                'init': {

                    #parameters for a 'defined'-type Layer
                    'defined': {
                        #raster to use for the Layer
                        'rast':                   stable,
                        #point coordinates
                        'pts':                    None,
                        #point values
                        'vals':                   None,
                        #interpolation method {None, 'linear', 'cubic',
                        #'nearest'}
                        'interp_method':          None,

                        }, # <END> 'defined'

                    }, # <END> 'init'

                }, # <END> layer num. 1


            #layer name (LAYER NAMES MUST BE UNIQUE!)
            'K': {

        #-------------------------------------#
        #--- layer num. 2: init parameters ---#
        #-------------------------------------#

                #initiating parameters for this layer
                'init': {

                    #parameters for a 'defined'-type Layer
                    'defined': {
                        #raster to use for the Layer
                        'rast':                   K,
                        #point coordinates
                        'pts':                    None,
                        #point values
                        'vals':                   None,
                        #interpolation method {None, 'linear', 'cubic',
                        #'nearest'}
                        'interp_method':          None,

                        }, # <END> 'defined'

                    }, # <END> 'init'

                }, # <END> layer num. 2


        #layer name (LAYER NAMES MUST BE UNIQUE!)
        'move': {

        #-------------------------------------#
        #--- layer num. 2: init parameters ---#
        #-------------------------------------#

                #initiating parameters for this layer
                'init': {

                    #parameters for a 'defined'-type Layer
                    'defined': {
                        #raster to use for the Layer
                        'rast':                   np.ones((50,50)),
                        #point coordinates
                        'pts':                    None,
                        #point values
                        'vals':                   None,
                        #interpolation method {None, 'linear', 'cubic',
                        #'nearest'}
                        'interp_method':          None,

                        }, # <END> 'defined'

                    }, # <END> 'init'

                }


    #### NOTE: Individual Layers' sections can be copy-and-pasted (and
    #### assigned distinct keys and names), to create additional Layers.


            } # <END> 'layers'

        }, # <END> 'landscape'


###############################################################################

###################
#### COMMUNITY ####
###################
    'comm': {

        'species': {

            #species name (SPECIES NAMES MUST BE UNIQUE!)
            'spp_0': {

            #-----------------------------------#
            #--- spp num. 0: init parameters ---#
            #-----------------------------------#

                'init': {
                    #starting number of individs
                    'N':                1000,
                    #carrying-capacity Layer name
                    'K_layer':          'K',
                    #multiplicative factor for carrying-capacity layer
                    'K_factor':         50,
                    }, # <END> 'init'

            #-------------------------------------#
            #--- spp num. 0: mating parameters ---#
            #-------------------------------------#

                'mating'    : {
                    #age(s) at sexual maturity (if tuple, female first)
                    'repro_age':                0,
                    #whether to assign sexes
                    'sex':                      False,
                    #ratio of males to females
                    'sex_ratio':                1/1,
                    #whether P(birth) should be weighted by parental dist
                    'dist_weighted_birth':       False,
                    #intrinsic growth rate
                    'R':                        0.5,
                    #intrinsic birth rate (MUST BE 0<=b<=1)
                    'b':                        0.5,
                    #expectation of distr of n offspring per mating pair
                    'n_births_distr_lambda':    1,
                    #whether n births should be fixed at n_births_dist_lambda
                    'n_births_fixed':           True,
                    #radius of mate-search area
                    'mating_radius':            5,
                    'choose_nearest_mate':      False,
                    'inverse_dist_mating':      False,
                    }, # <END> 'mating'

            #----------------------------------------#
            #--- spp num. 0: mortality parameters ---#
            #----------------------------------------#

                'mortality'     : {
                    #maximum age
                    'max_age':                      None,
                    #min P(death) (MUST BE 0<=d_min<=1)
                    'd_min':                        0,
                    #max P(death) (MUST BE 0<=d_max<=1)
                    'd_max':                        1,
                    #width of window used to estimate local pop density
                    'density_grid_window_width':    None,
                    }, # <END> 'mortality'

            #---------------------------------------#
            #--- spp num. 0: movement parameters ---#
            #---------------------------------------#

                'movement': {
                    #whether or not the species is mobile
                    'move':                     True,
                    #mode of distr of movement direction
                    'direction_distr_mu':       0,
                    #concentration of distr of movement direction
                    'direction_distr_kappa':    0,
                    #mean of distr of movement distance
                    'movement_distance_distr_param1':        0.25,
                    #variance of distr of movement distance
                    'movement_distance_distr_param2':     0.5,
                    'movement_distance_distr':             'wald',
                    #mean of distr of dispersal distance
                    'dispersal_distance_distr_param1':       0.5,
                    #variance of distr of dispersal distance
                    'dispersal_distance_distr_param2':    0.5,
                    'dispersal_distance_distr':             'wald',
                    #TODO: UNCOMMENT move_surf SECTION IF NEEDED!
                    #'move_surf':    {
                    #    'layer': 'move',
                    #    'mixture': True,
                    #    'vm_distr_kappa': 12,
                    #    'approx_len': 5000},
                    },    # <END> 'movement'


            #---------------------------------------------------#
            #--- spp num. 0: genomic architecture parameters ---#
            #---------------------------------------------------#

                'gen_arch': {
                    #file defining custom genomic arch
                    'gen_arch_file':            None,
                    #num of loci
                    'L':                        1000,
                    #value to use for fixed starting allele freqs (None to draw)
                    'start_p_fixed':            0.5,
                    #whether to start neutral locus freqs at 0
                    'start_neut_zero':          True,
                    #genome-wide per-base neutral mut rate (0 to disable)
                    'mu_neut':                  0,
                    #genome-wide per-base deleterious mut rate (0 to disable)
                    'mu_delet':                 0,
                    #shape of distr of deleterious effect sizes
                    'delet_alpha_distr_shape':  0.2,
                    #scale of distr of deleterious effect sizes
                    'delet_alpha_distr_scale':  0.2,
                    #alpha of distr of recomb rates
                    'r_distr_alpha':            1000,
                    #beta of distr of recomb rates
                    'r_distr_beta':             1e3,
                    #whether loci should be dominant (for allele '1')
                    'dom':                      False,
                    #whether to allow pleiotropy
                    'pleiotropy':               False,
                    #custom fn for drawing recomb rates
                    'recomb_rate_custom_fn':    None,
                    #number of recomb paths to hold in memory
                    'n_recomb_paths_mem':       int(1e4),
                    #total number of recomb paths to simulate
                    'n_recomb_paths_tot':       int(1e5),
                    'n_recomb_sims':            100_000,
                    #whether to generate recombination paths at each timestep
                    'allow_ad_hoc_recomb':       False,
                    #whether to save mutation logs
                    'mut_log':                  False,
                    'jitter_breakpoints': False,
                    'use_tskit': True,
                    'tskit_simp_interval': 100,

                    'traits': {

                        #-------------------------#
                        #---trait 0 parameters ---#
                        #-------------------------#
                        #trait name (TRAIT NAMES MUST BE UNIQUE!)
                        'trait_0': {
                            #trait-selection Layer name
                            'layer':                'shift',
                            #polygenic selection coefficient
                            'phi':                  1,
                            #number of loci underlying trait
                            'n_loci':               50,
                            #mutation rate at loci underlying trait
                            'mu':                   0,
                            #mean of distr of effect sizes
                            'alpha_distr_mu' :      0,
                            #variance of distr of effect size
                            'alpha_distr_sigma':    0,
                            #max allowed magnitude for an alpha value
                            'max_alpha_mag':        None,
                            #curvature of fitness function
                            'gamma':                1,
                            #whether the trait is universally advantageous
                            'univ_adv':             False
                            }, # <END> trait 0

                        #-------------------------#
                        #---trait 1 parameters ---#
                        #-------------------------#
                        #trait name (TRAIT NAMES MUST BE UNIQUE!)
                        'trait_1': {
                            #trait-selection Layer name
                            'layer':                'stable',
                            #polygenic selection coefficient
                            'phi':                  1,
                            #number of loci underlying trait
                            'n_loci':               50,
                            #mutation rate at loci underlying trait
                            'mu':                   0,
                            #mean of distr of effect sizes
                            'alpha_distr_mu' :      0,
                            #variance of distr of effect size
                            'alpha_distr_sigma':    0,
                            #max allowed magnitude for an alpha value
                            'max_alpha_mag':        None,
                            #curvature of fitness function
                            'gamma':                1,
                            #whether the trait is universally advantageous
                            'univ_adv':             False
                            }, # <END> trait 1


    #### NOTE: Individual Traits' sections can be copy-and-pasted (and
    #### assigned distinct keys and names), to create additional Traits.


                        }, # <END> 'traits'

                    }, # <END> 'gen_arch'


                }, # <END> spp num. 0



    #### NOTE: individual Species' sections can be copy-and-pasted (and
    #### assigned distinct keys and names), to create additional Species.


            }, # <END> 'species'

        }, # <END> 'comm'


###############################################################################

###############
#### MODEL ####
###############
    'model': {
        #total Model runtime (in timesteps)
        'T':            1250,
        #min burn-in runtime (in timesteps)
        'burn_T':       30,
        #seed number
        'num':          None,

        ###############################
        #### iterations parameters ####
        ###############################
        'its': {
            #num iterations
            'n_its':            1,
            #whether to randomize Landscape each iteration
            'rand_landscape':   False,
            #whether to randomize Community each iteration
            'rand_comm':        False,
            'rand_genarch':     True,
            #whether to burn in each iteration
            'repeat_burn':      False,
            }, # <END> 'iterations'


        ####################################
        #### data-collection parameters ####
        ####################################
        'data': {
            'sampling': {
                #sampling scheme {'all', 'random', 'point', 'transect'}
                'scheme':               'all',
                #sample size at each point, for point & transect sampling
                'n':                    1000,
                #coords of collection points, for point sampling
                'points':               None,
                #coords of transect endpoints, for transect sampling
                'transect_endpoints':   None,
                #num points along transect, for transect sampling
                'n_transect_points':    None,
                #collection radius around points, for point & transect sampling
                'radius':               None,
                #when to collect data
                'when':                 [1499, 1624, 1749],
                #whether to save current Layers when data is collected
                'include_landscape':    False,
                #whether to include fixed loci in VCF files
                'include_fixed_sites':  True,
                },
            'format': {
                #format for genetic data {'vcf', 'fasta'}
                'gen_format':           ['vcf'],
                #format for vector geodata {'csv', 'shapefile', 'geojson'}
                'geo_vect_format':      'csv',
                #format for raster geodata {'geotiff', 'txt'}
                'geo_rast_format':      'geotiff',
                'nonneut_loc_format':      'csv',
                },
            }, #<END> 'data'


        ######################################
        ##### stats-collection parameters ####
        ######################################
        #'stats': {
        #    #number of individs at time t
        #    'Nt': {
        #        #whether to calculate
        #        'calc':     False,
        #        #calculation frequency (in timesteps)
        #        'freq':     1,
        #        },
        #    #heterozgosity
        #    'het': {
        #        #whether to calculate
        #        'calc':     False,
        #        #calculation frequency (in timesteps)
        #        'freq':     5,
        #        #whether to mean across sampled individs
        #        'mean': False,
        #        },
        #    #minor allele freq
        #    'maf': {
        #        #whether to calculate
        #        'calc':     False,
        #        #calculation frequency (in timesteps)
        #        'freq':     5,
        #        },
        #    #mean fitness
        #    'mean_fit': {
        #        #whether to calculate
        #        'calc':     False,
        #        #calculation frequency (in timesteps)
        #        'freq':     5,
        #        },
        #    #linkage disequilibirum
        #    'ld': {
        #        #whether to calculate
        #        'calc':     False,
        #        #calculation frequency (in timesteps)
        #        'freq':     100,
        #        },
        #    }, # <END> 'stats'

        } # <END> 'model'

    } # <END> params

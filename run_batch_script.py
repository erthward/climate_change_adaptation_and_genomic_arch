def run_batch_script(i):
    print('STARTING BATCH NUMBER %i...\n' % i)
    import os
    #os.system("python batch_script.py")
    os.system(("python /global/home/users/drewhart/genarch_and_envchange/"
    	        "climate_change_adaptation_and_genomic_arch/batch_script.py"))
    return

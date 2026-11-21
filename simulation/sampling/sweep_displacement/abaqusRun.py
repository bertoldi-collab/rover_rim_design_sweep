# Import Tool
import os
import inspect
import math
import numpy as np
# Import tooling
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
samplingdir = os.path.dirname(currentdir)
simulationsdir = os.path.dirname(samplingdir)



#########################################################
#                       IMPORTS                         #
#########################################################
for root, dirs, files in os.walk('.'):
    for filename in files:
        if filename.endswith('.lck'):
            filepath = os.path.join(root, filename)
            os.remove(filepath)

# Import Abaqus classes
execfile(simulationsdir + '\\abaqus\\abaqus_rim_class.py')
# Import materials for wood and rubber
execfile(simulationsdir + '\\materials\\basswood_veneer.py')
matPar_basewood = matPar
execfile(simulationsdir + '\\materials\\rubber_band.py')
matPar_rubber = matPar

#########################################################
#                      PARAMETERS                       #
#########################################################
# Beam Stiffness values
hydration_E_table = {0: 11470.38, 1: 9165.85, 3: 5206.06, 7: 4326.80}

tasks = [
    (0.30, 7, True),
    (0.15, 7, True),
    (0.30, 0, True),
    (0.15, 0, True),
    (0.30, 1, True),
    (0.30, 3, True),
]

# Rubber stretch values
eps = (0.10, 0.20, 0.30, 0.40)

# Rubber thickness values
thicknesses = np.linspace(0.15, 0.50, 8)



# Rim parameters
radius = 62.5 # mm
beam_width = 5.0 # mm
mesh_seed = 0.490875  # approx 200 elements along beam length 
max_inward_actuation_mm = 25.0
arc_start_deg = 90.0 * 0.275 # starting angle for rubber band
arc_total_deg = 90.0 * 0.45  # ending angle for rubber band
theta_deg = arc_start_deg + arc_total_deg / 2.0
theta = math.radians(theta_deg)

# Simulation controller flags
overwrite_existing = False
reuse_existing_odb = True

#########################################################
#                         MODEL                         #
#########################################################
n_errors = 0
for wood_thickness, hydration_time, no_rubber_reference in tasks:
    hydration_E = hydration_E_table[hydration_time]
    hydration_label = 'hydration_{}min'.format(hydration_time)
    matPar_basewood['E'] = hydration_E
    resultsdir = os.path.join(
        currentdir, 'results', 'wood_{:.3f}mm'.format(wood_thickness))
    print('')
    print('CASE: wood {:.3f} mm, {} min hydration, E = {:.2f} MPa'.format(
        wood_thickness, hydration_time, hydration_E))


    #########################################################
    #                  WOOD-ONLY REFERENCE                  #
    #########################################################
    # Wood only reference case for future rubber comparison
    reference_dir = os.path.join(resultsdir, 'no_rubber', hydration_label)
    if no_rubber_reference and not os.path.exists(
            os.path.join(reference_dir, 'actuation.csv')):
        if not os.path.exists(reference_dir):
            os.makedirs(reference_dir)
        os.chdir(reference_dir)
        print('wood-only reference curve')

        Mdb()
        model = mdb.models['Model-1']
        name_relax = 'step-relax'
        name_actuation = 'step-actuation'
        model.StaticStep(name=name_relax, previous='Initial', nlgeom=ON,
            initialInc=0.001, minInc=1e-12, maxInc=0.01, maxNumInc=10000,
            timePeriod=1.0)

        sample = rim(model)
        sample.generate_geometry('wood_rim', radius)
        sample.material_properties(matPar_basewood, 'wood_rim', beam_width,
            wood_thickness, model.parts['wood_rim'].faces, matName='wood')
        sample.assembly('wood_rim_assem', 'wood_rim')
        sample.mesh('wood_rim', [B21], mesh_seed)
        sample.model.rootAssembly.regenerate()
        midpoint = model.rootAssembly.sets['wood_rim_assem.mid_point']
        sample.fixed_BC(None, name_relax, part_name='wood_rim_assem')

        # History output for relaxation step
        model.HistoryOutputRequest(name='Relax-Disp', createStepName=name_relax,
            variables=('U1', 'U2'), region=midpoint)
        model.HistoryOutputRequest(name='Relax-Force', createStepName=name_relax,
            variables=('RF1', 'RF2'), region=midpoint)
        model.HistoryOutputRequest(name='Relax-Energy', createStepName=name_relax,
            region=MODEL, variables=('ALLSE', 'ALLWK', 'ALLIE', 'ALLAE'))

        # Displacement-controlled Standard step, same as the rubber runs
        model.StaticStep(name=name_actuation, previous=name_relax, nlgeom=ON,
            initialInc=0.001, minInc=1e-10, maxInc=0.01, maxNumInc=50000,
            timePeriod=1.0, stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
            stabilizationMagnitude=0.0002)
        model.DisplacementBC(name='BC-RadialActuation',
            createStepName=name_actuation, region=midpoint,
            distributionType=UNIFORM, fieldName='', localCsys=None,
            u1=-max_inward_actuation_mm * math.cos(theta),
            u2=-max_inward_actuation_mm * math.sin(theta),
            u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
        model.HistoryOutputRequest(name='Actuation-Disp',
            createStepName=name_actuation, variables=('U1', 'U2'), region=midpoint)
        model.HistoryOutputRequest(name='Actuation-Force',
            createStepName=name_actuation, variables=('RF1', 'RF2'), region=midpoint)
        model.HistoryOutputRequest(name='Actuation-Energy',
            createStepName=name_actuation, region=MODEL,
            variables=('ALLSE', 'ALLWK', 'ALLIE', 'ALLAE'))

        sample.create_job('actuation', nCores=1)
        sim_time = sample.run_job()
        sample.post_process('actuation', name_actuation)
        os.chdir(currentdir)

    #########################################################
    #                        SWEEP                          #
    #########################################################
    variations = [(ep_a, t_rubber) for ep_a in eps for t_rubber in thicknesses]

    for ep_a, t_rubber in variations:
        case_dir = os.path.join(resultsdir,
            'rubber_thickness_{:.3f}mm'.format(t_rubber), hydration_label,
            'ep_a_{:.2f}'.format(ep_a))
        csv_file = os.path.join(case_dir, 'actuation.csv')
        odb_file = os.path.join(case_dir, 'Job.odb')
        failed_marker = os.path.join(case_dir, 'FAILED.txt')
        print('wood {:.3f} mm, {} min, rubber {:.3f} mm, ep_a {:.2f}'.format(
            wood_thickness, hydration_time, t_rubber, ep_a))

        # Skip what is already done, and re-read an ODB rather than re-solving
        if os.path.exists(csv_file) and not overwrite_existing:
            continue
        if not os.path.exists(case_dir):
            os.makedirs(case_dir)
        os.chdir(case_dir)
        if os.path.exists(odb_file) and reuse_existing_odb and not overwrite_existing:
            sample = rim(mdb.models['Model-1'])
            sample.jobName = 'Job'
            sample.post_process('actuation', 'step-accurate-actuation')
            os.chdir(currentdir)
            continue

        try:
            #########################################################
            #                       SIMULATION                      #
            #########################################################
            Mdb()
            model = mdb.models['Model-1']
            name_relax = 'step-relax'
            name_actuation = 'step-accurate-actuation'
            model.StaticStep(name=name_relax, previous='Initial', nlgeom=ON,
                initialInc=0.001, minInc=1e-12, maxInc=0.01, maxNumInc=10000,
                timePeriod=1.0)

            #########################################################
            #                       APPLICATION                     #
            #########################################################
            # Rubber band sits one half-thickness plus a 0.05 mm clearance
            # outside the wood ring
            gap = 0.5 * (wood_thickness + t_rubber) + 0.05

            sample = rim(model)
            sample.generate_geometry('wood_rim', radius)
            sample.generate_geometry('rubber_rim', radius + gap,
                arc_angle_total=arc_total_deg, arc_angle_start=arc_start_deg)
            sample.material_properties(matPar_basewood, 'wood_rim', beam_width,
                wood_thickness, model.parts['wood_rim'].faces, matName='wood')
            sample.material_properties(matPar_rubber, 'rubber_rim', beam_width,
                t_rubber, model.parts['rubber_rim'].faces, matName='rubber')
            sample.assembly('wood_rim_assem', 'wood_rim')
            sample.assembly('rubber_rim_assem', 'rubber_rim')
            sample.mesh('wood_rim', [B21], mesh_seed)
            sample.mesh('rubber_rim', [B21], mesh_seed)
            sample.model.rootAssembly.regenerate()
            midpoint = model.rootAssembly.sets['wood_rim_assem.mid_point']

            name_amp = 'amp-relax'
            sample.amplitude(name_amp, step_name=name_relax, num_steps=200,
                exp_fact=1, num_csv_TP=200, num_odb_TP=200, maxU=25)
            sample.fixed_BC(name_amp, name_relax, part_name='wood_rim_assem')
            sample.rim_EQ(name_amp, name_relax, part_name1='wood_rim_assem',
                part_name2='rubber_rim_assem')
            sample.rim_contact(name_relax, part_name1='wood_rim_assem',
                part_name2='rubber_rim_assem')
            # Stretched rubber band as a compressive eigenstrain
            sample.prestrain_BC(part_name='rubber_rim_assem', prestrain=ep_a,
                step=name_relax, amp=name_amp)
            model.HistoryOutputRequest(name='Relax-Disp', createStepName=name_relax,
                variables=('U1', 'U2'), region=midpoint)
            model.HistoryOutputRequest(name='Relax-Force', createStepName=name_relax,
                variables=('RF1', 'RF2'), region=midpoint)
            model.HistoryOutputRequest(name='Relax-Energy', createStepName=name_relax,
                region=MODEL, variables=('ALLSE', 'ALLWK', 'ALLIE', 'ALLAE'))

            # A single displacement-controlled Standard step preserves the stable
            # switched branch and keeps the whole actuation path in this ODB.
            model.StaticStep(name=name_actuation, previous=name_relax, nlgeom=ON,
                initialInc=0.001, minInc=1e-10, maxInc=0.01, maxNumInc=50000,
                timePeriod=1.0, stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
                stabilizationMagnitude=0.0002)
            model.DisplacementBC(name='BC-AccurateActuation',
                createStepName=name_actuation, region=midpoint,
                distributionType=UNIFORM, fieldName='', localCsys=None,
                u1=-max_inward_actuation_mm * math.cos(theta),
                u2=-max_inward_actuation_mm * math.sin(theta),
                u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
            model.HistoryOutputRequest(name='AccurateActuation-Disp',
                createStepName=name_actuation, variables=('U1', 'U2'), region=midpoint)
            model.HistoryOutputRequest(name='AccurateActuation-Force',
                createStepName=name_actuation, variables=('RF1', 'RF2'), region=midpoint)
            model.HistoryOutputRequest(name='AccurateActuation-Energy',
                createStepName=name_actuation, region=MODEL,
                variables=('ALLSE', 'ALLWK', 'ALLIE', 'ALLAE'))

            sample.create_job('Job', nCores=1)
            sim_time = sample.run_job()
            if not os.path.exists('Job.odb'):
                raise RuntimeError('Job ODB was not created.')
            sample.post_process('actuation', name_actuation)
            # A previous failure that has now succeeded should not keep its marker
            if os.path.exists(failed_marker):
                os.remove(failed_marker)
        except Exception as error:
            output = open(failed_marker, 'w')
            output.write(str(error) + '\n')
            output.close()
            n_errors = n_errors + 1
        os.chdir(currentdir)


print('')
print('Finished with {} errors'.format(n_errors))

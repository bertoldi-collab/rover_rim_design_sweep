from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import csv
import math
import numpy as np
import assembly 
import mesh
import time
import os
try:
    from odbAccess import openOdb
except:
    pass

#########################################################
#                   ABAQUS CLASS                        #
#########################################################

#   Basic abaqus class to generate a quarter rim object
#   for simulations in Abaqus 2025. 

class rim:
    def __init__(self,model):
        # Coordinate system for easier copy/paste from journal/replay files
        session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
        # Save model
        self.model = model


    def generate_geometry(self,name_part,R,arc_angle_total = 90.0,arc_angle_start=0.0):
        # Store relevant class parameters for rim
        self.arc_angle_total = arc_angle_total
        self.arc_angle_start = arc_angle_start
        self.R = R

        # Calculate rim parameters
        theta1 = math.radians(arc_angle_start)
        theta2 = math.radians(arc_angle_start + arc_angle_total)
        theta_mid = math.radians(arc_angle_start + arc_angle_total / 2.0)
        x1 = R * math.cos(theta1)
        y1 = R * math.sin(theta1)
        x2 = R * math.cos(theta2)
        y2 = R * math.sin(theta2)
        x_mid = R * math.cos(theta_mid)
        y_mid = R * math.sin(theta_mid)

        # Points for Rim object
        point1 = (x1, y1)
        point2 = (x2, y2)
        point_mid = (x_mid, y_mid)


        # Begin sketch
        self.model.ConstrainedSketch(name=name_part + '_sketch', sheetSize=200.0)
        # Create two arcs: point1 to midpoint, then midpoint to point2
        self.model.sketches[name_part + '_sketch'].ArcByCenterEnds(center=(0.0, 0.0)
            , direction=COUNTERCLOCKWISE, point1=point1, point2=point_mid)
        self.model.sketches[name_part + '_sketch'].ArcByCenterEnds(center=(0.0, 0.0)
            , direction=COUNTERCLOCKWISE, point1=point_mid, point2=point2)
        self.model.Part(dimensionality=TWO_D_PLANAR, name=name_part, type=
            DEFORMABLE_BODY)
        # Wire part
        self.model.parts[name_part].BaseWire(sketch=
            self.model.sketches[name_part + '_sketch'])
        del self.model.sketches[name_part + '_sketch']


        # Save vertices for future node reference
        endVertices = self.model.parts[name_part].vertices
        verts = (
            endVertices.findAt(coordinates=((point1[0], point1[1], 0.0),), ),
            endVertices.findAt(coordinates=((point2[0], point2[1], 0.0),), ))
        self.model.parts[name_part].Set(name='end_points', vertices=verts)

        # Save middle vertex for to find not for future loading
        midVertices = self.model.parts[name_part].vertices
        midVert = midVertices.findAt(coordinates=((point_mid[0], point_mid[1], 0.0),), )
        self.model.parts[name_part].Set(name='mid_point', vertices=(midVert,))



    def material_properties(self,matPar,name_part,w,t,faces,matName='Material'):
        # Abaqus Material Created
        self.model.Material(name=matName)

        # If Damping and Density provided, add them in
        try:
            self.model.materials[matName].Density(table=((matPar['rho'], ), ))
            self.model.materials[matName].Damping(alpha=matPar['alpha_damping'])
        except:
            print("Warning: Density and Damping parameters not implemented.")

        # Material property type
        if matPar['matModel'] == 'LIN':
            self.model.materials[matName].Elastic(table=((matPar['E'], matPar['nu']), ))

        if matPar['matModel'] == 'NH':
            self.model.materials[matName].Hyperelastic(materialType=ISOTROPIC, 
                testData=OFF, type=NEO_HOOKE, 
                table=((matPar['C10'], matPar['D']), ))

        
        # Create Material Section
        self.model.RectangularProfile(a=w, b=t, name='profile_'+matName)
        self.model.BeamSection(beamSectionOffset=(0.0, 0.0), 
            consistentMassMatrix=False, integration=DURING_ANALYSIS, material=
            matName, name='section_'+matName, poissonRatio=0.5, profile='profile_'+matName,
            temperatureVar=LINEAR)

        # Apply section to part
        self.model.parts[name_part].Set(edges=self.model.parts[name_part].edges, name='edges_'+matName)

        self.model.parts[name_part].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=
            self.model.parts[name_part].sets['edges_'+matName], sectionName=
            'section_'+matName, thicknessAssignment=FROM_SECTION)

        # Orientate beam section, use beam rendering tool in CAE to double check
        self.model.parts[name_part].assignBeamSectionOrientation(
            region=self.model.parts[name_part].sets['edges_'+matName],
            method=N1_COSINES, n1=(0.0, 0.0, -1.0))
        
        print('Successfully created and assigned material properties for '+matName+' material.')

        
    def mesh(self, name_part, elemCode, mesh_size_tile, elemLibrary=STANDARD):
        # Save mesh tile size for latter reference if necessasry
        self.mesh_size_tile = mesh_size_tile



        # Set type of element as bar element and seed
        self.model.parts[name_part].setElementType(elemTypes=(ElemType(
            elemCode=elemCode[0], elemLibrary=STANDARD), ), regions=(
            self.model.parts[name_part].edges, ))
        self.model.parts[name_part].seedPart(deviationFactor=0.1, 
            minSizeFactor=0.1, size=mesh_size_tile)

        # Generate Mesh
        self.model.parts[name_part].generateMesh(seedConstraintOverride=ON)


    def assembly(self,name_assem,name_part):
        # Generate Assembly
        self.name_assem = name_assem
        self.model.rootAssembly.Instance(dependent=ON, name=self.name_assem, part=self.model.parts[name_part])

    def amplitude(self,amp,num_steps,exp_fact,num_csv_TP,num_odb_TP,step_name='step-',start_point=0,end_point=1,amp_tuple=0,maxU=25,offset=0):
        # Exp_fact is numbers from 1 (linear ramp) to 100+ (exponential), allows for the tuning of load
        #step_name = step_name+self.name_part
        self.num_steps = num_steps
        self.exp_fact = exp_fact

        if amp_tuple != 0:
            self.model.TabularAmplitude(data=amp_tuple, name=
                amp, timeSpan=STEP)
        else:
            amp_data = self.exp_amp(num_steps,exp_fact,start_point,end_point,step_name,offset)
            self.model.TabularAmplitude(data=amp_data, name=
                amp, timeSpan=STEP)

        #Increase minimum number of attempts from 5 to 30
        if maxU == maxU:
            self.model.steps[step_name].control.setValues(allowPropagation=OFF, 
                resetDefaultValues=OFF, timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0, 
                12.0, 25, 6.0, 3.0, 50.0))
        self.num_csv_TP = num_csv_TP
        self.num_odb_TP = num_odb_TP
        amp_csv_TP = np.linspace(0,1,num_csv_TP)
        amp_odb_TP = np.linspace(0,1,num_odb_TP)

        time_csv_TP = num_csv_TP*[0]
        time_odb_TP = num_odb_TP*[0]
        if exp_fact == 1:
            for i,t in enumerate(np.linspace(0,1,num_csv_TP)):
                time_csv_TP[i] = ((float(t)),)
            for i,t in enumerate(np.linspace(0,1,num_odb_TP)):
                time_odb_TP[i] = ((float(t)),) 
        else:
            for i in range(num_csv_TP):
                    if i != 0:
                        time_csv_TP[i] = (((np.log10(amp_csv_TP[i]*((exp_fact)-1)+1))/np.log10(exp_fact)),)
                    else:
                        time_csv_TP[i] = ((0),)
            for i in range(num_odb_TP):
                    if i != 0:
                        time_odb_TP[i] = (((np.log10(amp_odb_TP[i]*((exp_fact)-1)+1))/np.log10(exp_fact)),)
                    else:
                        time_odb_TP[i] = ((0),)
        self.model.TimePoint(name='TimePoints-CSV', points=tuple(time_csv_TP))
        self.model.TimePoint(name='TimePoints-ODB', points=tuple(time_odb_TP))

    def exp_amp(self,num_timesteps,b,start_point,end_point,step_name,offset):
        # Exponential amplitude function generated by Helen Read
        def exp_fxn(x,b):
            if b ==1:
                return x
            else:
                return 1./(b-1) * (b**x - 1)

        all_t = np.linspace(start_point,end_point,num_timesteps)
        all_tuples = (num_timesteps)*[0]
        for i in range(num_timesteps):
            all_tuples[i] = (all_t[i],exp_fxn(all_t[i],b)+offset)
        return all_tuples
    
    def fixed_BC(self, amp, step,part_name):
        A = self.model.rootAssembly

        # Fixed BC on either side of the rim
        self.model.DisplacementBC(
            createStepName=step, distributionType=UNIFORM,
            fieldName='', localCsys=None, name='BC-ENCASTRE',
            region=A.sets[part_name + '.end_points'],
            u1=0, u2=0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=0)

    def rim_EQ(self, amp, step, part_name1, part_name2):
        A = self.model.rootAssembly
        # Get coordinates of rubber endpoints
        rubber_top_node = None
        rubber_bottom_node = None
        print(A)
        print(part_name1, part_name2)
        print(A.instances[part_name1].nodes)
        for nodes in A.instances[part_name2].nodes:
            print(nodes.coordinates)
            if rubber_top_node is None or nodes.coordinates[1] > rubber_top_node.coordinates[1]:
                rubber_top_node = nodes

            if rubber_bottom_node is None or nodes.coordinates[1] < rubber_bottom_node.coordinates[1]:
                rubber_bottom_node = nodes
        # Get wood nodes that are closest to endpoints of rubber.
        wood_top_node = None
        wood_bottom_node = None
        min_dist_top = float('inf')
        min_dist_bottom = float('inf')

        # Minimize to connect the two closests nodes
        for node in A.instances[part_name1].nodes:
            dist_to_top = np.linalg.norm(np.array(node.coordinates) - np.array(rubber_top_node.coordinates))
            dist_to_bottom = np.linalg.norm(np.array(node.coordinates) - np.array(rubber_bottom_node.coordinates))
            
            if dist_to_top < min_dist_top:
                wood_top_node = node
                min_dist_top = dist_to_top
            
            if dist_to_bottom < min_dist_bottom:
                wood_bottom_node = node
                min_dist_bottom = dist_to_bottom


        # Eq sets
        A.Set(name='wood_top_node', nodes=A.instances[part_name1].nodes.sequenceFromLabels((wood_top_node.label,)))
        A.Set(name='wood_bottom_node', nodes=A.instances[part_name1].nodes.sequenceFromLabels((wood_bottom_node.label,)))
        A.Set(name='rubber_top_node', nodes=A.instances[part_name2].nodes.sequenceFromLabels((rubber_top_node.label,)))
        A.Set(name='rubber_bottom_node', nodes=A.instances[part_name2].nodes.sequenceFromLabels((rubber_bottom_node.label,)))

        # Use equation to tie nodes together
        print(wood_top_node, rubber_top_node)
        for degF in [1, 2,6]:
            self.model.Equation(name='EQ-Top-{}'.format(degF), terms=((1.0, 'wood_top_node', degF, -1), (-1.0, 'rubber_top_node', degF, -1)))
  
        # Create equation constraint for bottom nodes
        for degF in [1, 2,6]:
            self.model.Equation(name='EQ-Bottom-{}'.format(degF), terms=((1.0, 'wood_bottom_node', degF, -1), (-1.0, 'rubber_bottom_node', degF, -1)))

    def rim_contact(self, step, part_name1, part_name2, friction=0.0):
        """Create element-based Abaqus/Standard contact for the two B21 rims."""
        A = self.model.rootAssembly
        wood = A.instances[part_name1]
        rubber = A.instances[part_name2]


        if len(wood.elements) == 0 or len(rubber.elements) == 0:
            raise RuntimeError('Mesh both contact instances before creating B21 contact.')

        wood_types = set(str(element.type).upper() for element in wood.elements)
        rubber_types = set(str(element.type).upper() for element in rubber.elements)
        if any(not element_type.startswith('B21') for element_type in wood_types):
            raise RuntimeError(
                '{} must contain only B21 elements; found {}.'
                .format(part_name1, sorted(wood_types)))
        if any(not element_type.startswith('B21') for element_type in rubber_types):
            raise RuntimeError(
                '{} must contain only B21 elements; found {}.'
                .format(part_name2, sorted(rubber_types)))

        def side1_faces_inward(instance):
            """Verify this model's B21 SIDE1 normal is radially inward."""
            nodes_by_label = dict(
                (node.label, node) for node in instance.nodes)
            for element in instance.elements:
                def resolve_node(connectivity_id):
                    if 0 <= connectivity_id < len(instance.nodes):
                        return instance.nodes[connectivity_id]
                    return nodes_by_label.get(connectivity_id)

                node1 = resolve_node(element.connectivity[0])
                node2 = resolve_node(element.connectivity[1])
                if node1 is None or node2 is None:
                    raise RuntimeError(
                        'Cannot resolve B21 connectivity {} in {}.'
                        .format(element.connectivity, instance.name))
                dx = node2.coordinates[0] - node1.coordinates[0]
                dy = node2.coordinates[1] - node1.coordinates[1]
                length = math.sqrt(dx * dx + dy * dy)
                if length == 0.0:
                    raise RuntimeError('Zero-length B21 element in {}.'.format(instance.name))

                # CAE displays SIDE1 as the in-plane normal -tangent x +Z.
                normal_x = -dy / length
                normal_y = dx / length
                mid_x = 0.5 * (node1.coordinates[0] + node2.coordinates[0])
                mid_y = 0.5 * (node1.coordinates[1] + node2.coordinates[1])
                if normal_x * mid_x + normal_y * mid_y >= 0.0:
                    return False
            return True


        # Existing arcs are concentric and counter-clockwise. 
        if not side1_faces_inward(wood):
            raise RuntimeError(
                'Cannot determine an inward SIDE1 for {} from B21 connectivity.'
                .format(part_name1))
        if not side1_faces_inward(rubber):
            raise RuntimeError(
                'Cannot determine an inward SIDE1 for {} from B21 connectivity.'
                .format(part_name2))

        for name in ('BeamBeamContact', 'GC-WoodRubber', 'Wood-Rubber-Contact'):
            if name in self.model.interactions.keys():
                del self.model.interactions[name]
        for name in ('BeamContactWood', 'BeamContactRubber', 'surf_wood', 'surf_rubber'):
            if name in A.surfaces.keys():
                del A.surfaces[name]
        for name in ('BeamContactProp', 'ContactProp-WoodRubber'):
            if name in self.model.interactionProperties.keys():
                del self.model.interactionProperties[name]

        # These are assembly MeshElementArray objects, not geometry Edge objects.
        wood_surface = A.Surface(
            name='BeamContactWood', side2Elements=wood.elements)
        rubber_surface = A.Surface(
            name='BeamContactRubber', side1Elements=rubber.elements)

        self.model.ContactProperty('BeamContactProp')
        self.model.interactionProperties['BeamContactProp'].NormalBehavior(
            pressureOverclosure=HARD,
            allowSeparation=ON,
            constraintEnforcementMethod=DEFAULT)


        # Wood is the stiffer main surface; rubber is the softer secondary.
        self.model.SurfaceToSurfaceContactStd(
            name='BeamBeamContact',
            createStepName='Initial',
            main=wood_surface,
            secondary=rubber_surface,
            sliding=FINITE,
            interactionProperty='BeamContactProp',
            enforcement=SURFACE_TO_SURFACE,
            thickness=ON)

    def prestrain_BC(self, part_name, prestrain, matName='rubber', step=None, amp=None):
        A = self.model.rootAssembly
        ep = float(prestrain)
        # Because we start from a stretched state, the compressive eigenstrain is negative.
        eps_c = -ep / (1.0 + ep)
        region = A.sets[part_name + '.edges_' + matName]



        self.model.materials[matName].Expansion(table=((1.0,),), zero=0.0)

        self.model.Temperature(
            name='Prestrain-Initial',
            createStepName='Initial',
            region=region,
            distributionType=UNIFORM,
            crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
            magnitudes=(0.0,))

        self.model.Temperature(
            name='Prestrain-Relax',
            createStepName=step,
            region=region,
            distributionType=UNIFORM,
            crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
            magnitudes=(eps_c,),
            amplitude=amp)
        return

    def create_time_points(self, num_csv_TP=50, num_odb_TP=None, time_period=1.0):
        # Time points for output settings
        if num_odb_TP is None:
            num_odb_TP = num_csv_TP
        self.num_csv_TP = int(num_csv_TP)
        self.num_odb_TP = int(num_odb_TP)
        self.output_time_period = float(time_period)

        for name in ('TimePoints-CSV', 'TimePoints-ODB'):
            if name in self.model.timePoints.keys():
                del self.model.timePoints[name]

        csv_pts = tuple(
            (float(t),) for t in np.linspace(0.0, float(time_period), int(num_csv_TP)))
        odb_pts = tuple(
            (float(t),) for t in np.linspace(0.0, float(time_period), int(num_odb_TP)))
        self.model.TimePoint(name='TimePoints-CSV', points=csv_pts)
        self.model.TimePoint(name='TimePoints-ODB', points=odb_pts)
        print('Created TimePoints-CSV ({}), TimePoints-ODB ({}).'.format(
            num_csv_TP, num_odb_TP))

    def create_job(self,jobName,nCores,file_subroutine=False,path_subroutine='',initial_conditions=False,Gshear=0,MPC=False,Corr_Normal=False):
        # Regenerate assembly
        self.model.rootAssembly.regenerate()
        # Patch MAXRATIO after regenerate so the keyword block stays current
        if hasattr(self, '_contact_surface_name'):
            self.add_surface_max_ratio(
                self._contact_surface_name,
                maxRatio=getattr(self, '_contact_max_ratio', 0.05))
        self.jobName = jobName
        self.file_subroutine = file_subroutine
        self.nCores = nCores
        # Create temporary job
        jobName_tmp = 'job_tmp'
        if self.file_subroutine == False:
            mdb.Job(
                name=jobName_tmp,
                numCpus=self.nCores, 
                numDomains=self.nCores, 
                explicitPrecision=DOUBLE_PLUS_PACK,
                nodalOutputPrecision=SINGLE, 
                atTime=None, 
                contactPrint=OFF, 
                description='', 
                echoPrint=OFF, 
                getMemoryFromAnalysis=True, 
                historyPrint=OFF, 
                memory=90, 
                memoryUnits=PERCENTAGE, 
                model='Model-1', 
                modelPrint=OFF, 
                multiprocessingMode=THREADS, 
                numGPUs=0, 
                parallelizationMethodExplicit=DOMAIN, 
                activateLoadBalancing=False, 
                queue=None, 
                resultsFormat=ODB, 
                scratch='', 
                type=ANALYSIS, 
                waitHours=0, 
                waitMinutes=0)

            # Update input file to include material orientations
            mdb.jobs[jobName_tmp].writeInput()

            oldFile = open(jobName_tmp + '.inp', 'r')
            newFile = open(self.jobName + '.inp', 'w+')
            
            # linesToSkip = 0
            for line in oldFile:
                # if linesToSkip > 0:
                #     linesToSkip -= 1
                #     continue

                if line.startswith('*Material,') and Gshear!=0:
                    newFile.write(line)
                    newFile.write('*Transverse Shear\n{},\n'.format(Gshear))
                elif line.startswith('0, MPC-corner') and MPC:
                    line = line[0:23]+'MPC-origin,'+line[23:-1]
                    newFile.write(line)
                else:
                    newFile.write(line)

            oldFile.close()
            newFile.close()

            # Create final job
            mdb.JobFromInputFile(
                name=self.jobName, 
                inputFileName=self.jobName+'.inp',
                numCpus=self.nCores,
                numDomains=self.nCores,
                explicitPrecision=DOUBLE_PLUS_PACK,
                nodalOutputPrecision=SINGLE, 
                atTime=None, 
                getMemoryFromAnalysis=True, 
                memory=90, 
                memoryUnits=PERCENTAGE, 
                multiprocessingMode=THREADS, 
                numGPUs=0, 
                parallelizationMethodExplicit=DOMAIN, 
                activateLoadBalancing=False, 
                queue=None, 
                resultsFormat=ODB, 
                scratch='', 
                type=ANALYSIS, 
                waitHours=0, 
                waitMinutes=0)
        else:
            mdb.Job(
                name=jobName_tmp,
                numCpus=self.nCores, 
                numDomains=self.nCores, 
                explicitPrecision=DOUBLE_PLUS_PACK,
                nodalOutputPrecision=SINGLE, 
                userSubroutine=path_subroutine,
                atTime=None, 
                contactPrint=OFF, 
                description='', 
                echoPrint=OFF, 
                getMemoryFromAnalysis=True, 
                historyPrint=OFF, 
                memory=90, 
                memoryUnits=PERCENTAGE, 
                model='Model-1', 
                modelPrint=OFF, 
                multiprocessingMode=THREADS, 
                numGPUs=0, 
                parallelizationMethodExplicit=DOMAIN, 
                activateLoadBalancing=False, 
                queue=None, 
                resultsFormat=ODB, 
                scratch='', 
                type=ANALYSIS, 
                waitHours=0, 
                waitMinutes=0)
            # Update input file to include material orientations
            mdb.jobs[jobName_tmp].writeInput()
            oldFile = open(jobName_tmp + '.inp', 'r')
            newFile = open(self.jobName + '.inp', 'w+')
            node_normal = False
            # linesToSkip = 0
            for line in oldFile:
                # if linesToSkip > 0:
                #     linesToSkip -= 1
                if node_normal:
                    if line.startswith('*Element'):
                        node_normal=False
                    else:
                        # try:
                        
                        # print(line_temp)
                        try:
                            line_temp = np.array(line.replace(',','').split(),dtype=float).tolist()
                            normal = np.array([line_temp[1],line_temp[2],0])/np.linalg.norm(np.array([line_temp[1],line_temp[2],0]))
                            # Convert NumPy array to list and add to numbers_list
                            line_temp.extend(normal.tolist())
                            line_temp[0] = int(line_temp[0])

                            # Convert the list back to a string with commas
                            line = ','.join(map(str, line_temp))
                            line = line+'\n'
 
                        except:
                            print('Zero normal encountered')
                            node_normal = False

                
                if line.startswith('*Orientation,'):
                    newFile.write(line)
                    # newFile.write(line[:-1] + ',local directions=2\n')
                elif line.startswith('*Transverse Shear'):
                    newFile.write('*Transverse Shear Stiffness\n')
                elif line =='*Node\n' and Corr_Normal:
                    node_normal = True
                    newFile.write(line)
                elif line.startswith('0, MPC-corner') and MPC:
                    line = line[0:22]+'MPC-origin,'+line[22:-1]+'\n'
                    newFile.write(line)
                elif line.startswith('0, MPC-node') and MPC:
                    line_temp = line.split(',')
                    node_temp = line_temp[1].split('-')
                    node_temp[2] = '1'
                    line_temp[-1] = '-'.join(node_temp)
                    line_temp.append('RefPoInt-Fluid-Cavity')
                    line = ','.join(line_temp)+'\n'
                    newFile.write(line)
                else:
                    newFile.write(line)

            oldFile.close()
            newFile.close()

            # Create final job
            mdb.JobFromInputFile(
                name=self.jobName, 
                inputFileName=self.jobName+'.inp',
                numCpus=self.nCores,
                numDomains=self.nCores,
                explicitPrecision=DOUBLE_PLUS_PACK,
                nodalOutputPrecision=SINGLE, 
                userSubroutine=path_subroutine,
                atTime=None, 
                getMemoryFromAnalysis=True, 
                memory=90, 
                memoryUnits=PERCENTAGE, 
                multiprocessingMode=THREADS, 
                numGPUs=0, 
                parallelizationMethodExplicit=DOMAIN, 
                activateLoadBalancing=False, 
                queue=None, 
                resultsFormat=ODB, 
                scratch='', 
                type=ANALYSIS, 
                waitHours=0, 
                waitMinutes=0)
        del mdb.jobs['job_tmp']

    def run_job(self, waitForCompletion=True, timeOut=6000, retainFiles=None):
        t1 = time.time()
        mdb.jobs[self.jobName].submit(consistencyChecking=OFF)
        print('Submitted: {}'.format(self.jobName))
        if waitForCompletion:
            try:
                print('Waiting: {}'.format(self.jobName))
                mdb.jobs[self.jobName].waitForCompletion(timeOut)
            except:
                print('Error: Job timed out or failed to complete.')
        t2 = time.time()
        t_total = t2 - t1
        time.sleep(2)
        if retainFiles is not None:
            self.clean_up(retainFiles=retainFiles)
        return t_total

    def clean_up(self, retainFiles=None):
        if retainFiles is None:
            retainFiles = ['.odb', '.sta', '.msg', '.inp']
        extensions = ['.odb', '.sta', '.dat', '.com', '.ipm', '.log', '.prt',
                      '.sim', '.msg', '.lck', '.inp', '.env', '.cid', '.pmg',
                      '.abq', '.pac', '.sel', '.stt', '.mdl', '.res', '.nsa',
                      '.psf', '.fil']
        for e in retainFiles:
            if e in extensions:
                extensions.remove(e)
        for e in extensions:
            try:
                os.remove(self.jobName + e)
            except:
                pass
        try:
            os.remove('job_tmp.inp')
        except:
            pass

    def post_process(self, fileName, step_name):
        # Save the tracked node history for one step: time, U1, U2, RF1, RF2
        odb = openOdb(self.jobName + '.odb')
        step = odb.steps[step_name]

        def history(variable):
            for name in step.historyRegions.keys():
                outputs = step.historyRegions[name].historyOutputs
                if variable in outputs.keys():
                    return outputs[variable].data
            raise RuntimeError('{} is missing from {}.'.format(variable, step_name))
        u1 = history('U1')
        u2 = history('U2')
        rf1 = history('RF1')
        rf2 = history('RF2')

        RES = np.zeros((len(u1), 5))
        for i in range(len(u1)):
            RES[i][0] = u1[i][0]
            RES[i][1] = u1[i][1]
            RES[i][2] = u2[i][1]
            RES[i][3] = rf1[i][1]
            RES[i][4] = rf2[i][1]

        np.savetxt(fileName + '.csv', RES, fmt='%.12g', delimiter=',',
            header='Time,U1,U2,RF1,RF2', comments='')
        odb.close()
        print('Post-processing is done. File saved: ' + fileName + '.csv')
        return

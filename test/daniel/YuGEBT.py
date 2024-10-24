# -*- coding: utf-8 -*-
# @Author: damov
# @Date:   2018-07-11T15:43:32+01:00
# @Last modified by:   damov
# @Last modified time: 2018-08-16T21:09:14+01:00

import os, re
import numpy as np
from fem_driver import FEMDriver #............................................. FEM driver abstract class

# Import FEM beam mesh
from ..FEM.Mesh import FEMesh
from ..FEM.ElementList import FEElementList
from ..FEM.NodeList import FENodeList
from falcons.include.cScreenOutput import ScreenOutput

"""
--------------------------------------------------------------------------------
    Class representing the FEM beam GEBT driver, implimented by Wenbin Yu1 and
    Qi Wang (https://cdmhub.org/resources/gebt)
--------------------------------------------------------------------------------
"""
class YuGEBT(FEMDriver):
    run_mode = ["steady", "unsteady"] #........................................ Modes in which the gebt solver can run


    def __init__(self, param, fem_mesh, message=ScreenOutput()):
        super(YuGEBT, self).__init__(param, fem_mesh, message=message) #....... Initialize super class
    #---Define axis orientation-------------------------------------------------
        self.axis_orientation = np.ones(6)
        self.axis_orientation[0] = 1.0
        self.axis_orientation[1] = 1.0
        self.axis_orientation[2] = 1.0
        self.axis_orientation[3] = 1.0
        self.axis_orientation[4] = 1.0
        self.axis_orientation[5] = 1.0
    #---Return------------------------------------------------------------------
        return

    """
    ----------------------------------------------------------------------------
        Initialize map
    ----------------------------------------------------------------------------
    """
    def Init(self):
        self.time_step = 0 #.................................................... Physical time step
        self.iteration = 0 #.................................................... Number of execution of the GEBT code

        self.analysis_type = self.param.Get("Solver mode") #.................... Analysis type
        self.inpt_filepref = self.param.Get("Prefix of input file of driver") #. Path where the input file has to be written
        self.niter         = self.param.Get("Number of gebt nonlinear iteration") #. Number of steps in the fem code
        self.nsteps        = self.param.Get("Number of gebt nonlinear steps")
        self.solpath_pref  = self.param.Get("Tecplot beam result output file prefix")

        self.delay_iter = self.param.Get_or_Default("Delay coupling by iterations", default=0)
        self.ramp_step  = self.param.Get_or_Default("Initial amount of time steps for force ramp gain", default=None)

        self.dt = 0.0
        if self.analysis_type == "unsteady":
            self.dt = self.param.Get("Unsteady physical time step size")

        self.relax_factor = self.param.Get("Structural displacement relaxation factor")

    #---Stop if mode type is unknown--------------------------------------------
        if not self.analysis_type in YuGEBT.run_mode:
            self.message.error("GEBT FEM can not run in following mode: %s" %(self.analysis_type))

    #---Matrices of initial values----------------------------------------------
        self.initial_elem_pos     = np.zeros((len(self.fem_mesh.elements),3)) #. Initial element position
        self.initial_elem_rot     = np.zeros((len(self.fem_mesh.elements),3)) #. Initial element rotation
        self.initial_elem_pos_vel = np.zeros((len(self.fem_mesh.elements),3)) #. Initial element linear velocity
        self.initial_elem_rot_vel = np.zeros((len(self.fem_mesh.elements),3)) #. Initial element angular velocity

    #---Initial conditions file-------------------------------------------------
        self.ini_conditions = np.zeros( (2*len(self.fem_mesh.elements),6) )

    #---Return------------------------------------------------------------------
        return

    """
    ----------------------------------------------------------------------------
        Run one steady-state iteration
    ----------------------------------------------------------------------------
    """
    def Run(self):
        self.iteration += 1 #........................................................ Increment iteration

    #---Update file path--------------------------------------------------------
        if self.analysis_type == "steady":
            self.file_path_input  = self.inpt_filepref + "_iter=" + str(self.iteration) #. file path to input file
        elif self.analysis_type == "unsteady":
            self.file_path_input  = self.inpt_filepref \
                                  + "_time_step=" + str(self.time_step) \
                                  + "_iter=" + str(self.iteration) #. file path to input file
        self.file_path_output = self.file_path_input + ".out" #....................... file path to output file

    #---Write input file--------------------------------------------------------
        self.__write_input_file()
        if self.analysis_type == "unsteady":
            self.__write_ini_file()

    #---Execute simulation------------------------------------------------------
        cmd = "gebt " + self.file_path_input #................................. Command string
        out = os.popen(cmd).read() #........................................... Execute command

    #---Check if solution converged---------------------------------------------
        error_msg = "The solution does not converge after the maximum number of iterations"
        if error_msg in out:
            raise Exception("GEBT: The solution does not converge after the maximum number of iterations")

    #---Check if there are Error messages---------------------------------------
        if "ERROR RETURN" in out:
            raise Exception("GEBT: There are error message in GEBT input file: %s " %(self.file_path_input))

    #---Print output------------------------------------------------------------
        num_gebt_iter = len(re.findall("ITERATION=.*\n", out))
        self.message.out("\n")
        self.message.out("Calculation of gebt finished in %i iterations\n" %(num_gebt_iter))
        self.message.out("\n")

    #---Read output file--------------------------------------------------------
        self.__read_output_file()

    #---return------------------------------------------------------------------
        return

    """
    ----------------------------------------------------------------------------
        Save initial state
    ----------------------------------------------------------------------------
    """
    def SaveInitialState(self):
        path_file_ini_out = self.file_path_input + ".ini_out"
        with open(path_file_ini_out,"r") as f:
            for k, line in enumerate(f):
                self.ini_conditions[k,:] = np.array([float(v) for v in line.split() if v])
                #if k >= len(self.fem_mesh.elements):
                #    self.ini_conditions[k,:] = np.array([float(v) for v in line.split() if v])
                #else:
                #    self.ini_conditions[k,:] = np.zeros(6)
    #---Return------------------------------------------------------------------
        return


    """
    ----------------------------------------------------------------------------
        Function which performs after subiterations
    ----------------------------------------------------------------------------
    """
    """
    def AfterUnsteadyIter(self):
        self.time_step += 1
        self.iteration = 0
        self.__read_inertial_values() #........................................ Read linear and angular element velocities
    #---Return------------------------------------------------------------------
        return
    """

    """
    ----------------------------------------------------------------------------
        Dump solution of beam
    ----------------------------------------------------------------------------
    """
    def Dump_solution(self, title=None):
    #---Create file path--------------------------------------------------------
        if self.analysis_type == "steady":
            path = self.solpath_pref + "_iter=" + str(self.iteration)
            time = self.iteration
        elif self.analysis_type == "unsteady":
            path = self.solpath_pref + "_time_step=" + str(self.time_step)
            time = self.time_step
    #---Dump fem mesh to solution-----------------------------------------------
        self.fem_mesh.DumpTecplot(path, title=title, time=time)
    #---Return------------------------------------------------------------------
        return

    """
    ----------------------------------------------------------------------------
        Read output file
    ----------------------------------------------------------------------------
    """
    def __read_output_file(self):
    #---Read output file--------------------------------------------------------
        f = open(self.file_path_output, "r")
        output_txt = f.read()
        f.close()

    #---Reset all velocities----------------------------------------------------
        for node in self.fem_mesh.nodes:
            node.v1     = 0.0 #................................................ Reset linear velocity of the element at mid point in global coordinate system a
            node.v2     = 0.0 #................................................ Reset linear velocity of the element at mid point in global coordinate system a
            node.v3     = 0.0 #................................................ Reset linear velocity of the element at mid point in global coordinate system a
            node.theta1 = 0.0 #................................................ Reset angular velocity of the element at mid point in global coordinate system a
            node.theta2 = 0.0 #................................................ Reset angular velocity of the element at mid point in global coordinate system a
            node.theta3 = 0.0 #................................................ Reset angular velocity of the element at mid point in global coordinate system a

    #---Read point values and save to nodes-------------------------------------
        node_data = re.findall(r"Point.*\n.*\n.*\n.*\n.*\n", output_txt) #..... Get from output file all point data
        for i, nd_str in enumerate(node_data):
            node_id = int(re.findall(r"Point #:\s*\d+", nd_str)[0].split()[-1])
            values_nd = [float(v) for v in re.findall(r"-?\d*\.\d*[eE][+-]\d*", nd_str)]

            #print(node_id, len(values_nd))

        #---Unpack values-------------------------------------------------------
            x0 =                            values_nd[0] #...................... Undeformed position of the node in x-direction
            y0 =                            values_nd[1] #...................... Undeformed position of the node in y-direction
            z0 =                            values_nd[2] #...................... Undeformed position of the node in z-direction
            u1 = self.axis_orientation[0] * values_nd[3] #...................... Displacement of the node in x-direction
            u2 = self.axis_orientation[1] * values_nd[4] #...................... Displacement of the node in y-direction
            u3 = self.axis_orientation[2] * values_nd[5] #...................... Displacement of the node in z-direction
            t1 = self.axis_orientation[3] * values_nd[6] #...................... Rotation around x-axis
            t2 = self.axis_orientation[4] * values_nd[7] #...................... Rotation around y-axis
            t3 = self.axis_orientation[5] * values_nd[8] #...................... Rotation around z-axis
            F1 = self.axis_orientation[0] * values_nd[9] #...................... Inner force in x-direction
            F2 = self.axis_orientation[1] * values_nd[10] #..................... Inner force in y-direction
            F3 = self.axis_orientation[2] * values_nd[11] #..................... Inner force in z-direction
            M1 = self.axis_orientation[3] * values_nd[12] #..................... Inner moment around x-axis
            M2 = self.axis_orientation[4] * values_nd[13] #..................... Inner moment around y-axis
            M3 = self.axis_orientation[5] * values_nd[14] #..................... Inner moment around z-axis

        #---Pick the node by coordinate-----------------------------------------
            #node = self.fem_mesh.nodes.GetNearestToXYZ(x0, y0, z0)
            node = self.fem_mesh.nodes[node_id-1]

            node.dx    = node.dx    + self.relax_factor * (u1 - node.dx   ) #... Displacement in x direction
            node.dy    = node.dy    + self.relax_factor * (u2 - node.dy   ) #... Displacement in y direction
            node.dz    = node.dz    + self.relax_factor * (u3 - node.dz   ) #... Displacement in z direction
            node.dphix = node.dphix + self.relax_factor * (t1 - node.dphix) #... Twist incriment of the node in x direction
            node.dphiy = node.dphiy + self.relax_factor * (t2 - node.dphiy) #... Twist incriment of the node in y direction
            node.dphiz = node.dphiz + self.relax_factor * (t3 - node.dphiz) #... Twist incriment of the node in z direction

    #---Read data at mid points of of elements and save internal beam values ---
        if self.analysis_type == "unsteady":
            elem_data = re.findall(r"Member.*\n.*\n.*\n.*\n.*\n.*\n", output_txt)
        else:
            elem_data = re.findall(r"Member.*\n.*\n.*\n.*\n.*\n", output_txt)

        for i, nd_str in enumerate(elem_data):
            member_id = int(re.findall(r"Member.*\n", nd_str)[0].split()[-1])
            values_nd = [float(v) for v in re.findall(r"-?\d*\.\d*[eE][+-]\d*", nd_str)]

            x0 =                            values_nd[0] #..................... Undeformed position of the element midpoint in x-direction
            y0 =                            values_nd[1] #..................... Undeformed position of the element midpoint in y-direction
            z0 =                            values_nd[2] #..................... Undeformed position of the element midpoint in z-direction
            u1 = self.axis_orientation[0] * values_nd[3] #..................... Displacement of the node in x-direction
            u2 = self.axis_orientation[1] * values_nd[4] #..................... Displacement of the node in y-direction
            u3 = self.axis_orientation[2] * values_nd[5] #..................... Displacement of the node in z-direction
            t1 = self.axis_orientation[3] * values_nd[6] #..................... Rotation around x-axis
            t2 = self.axis_orientation[4] * values_nd[7] #..................... Rotation around y-axis
            t3 = self.axis_orientation[5] * values_nd[8] #..................... Rotation around z-axis
            F1 = self.axis_orientation[0] * values_nd[9] #..................... Inner force in x-direction
            F2 = self.axis_orientation[1] * values_nd[10] #.................... Inner force in y-direction
            F3 = self.axis_orientation[2] * values_nd[11] #.................... Inner force in z-direction
            M1 = self.axis_orientation[3] * values_nd[12] #.................... Inner moment around x-axis
            M2 = self.axis_orientation[4] * values_nd[13] #.................... Inner moment around y-axis
            M3 = self.axis_orientation[5] * values_nd[14] #.................... Inner moment around z-axis

            if self.analysis_type == "unsteady":
                P1 = self.axis_orientation[0] * values_nd[15]
                P2 = self.axis_orientation[1] * values_nd[16]
                P3 = self.axis_orientation[2] * values_nd[17]
                H1 = self.axis_orientation[3] * values_nd[18]
                H2 = self.axis_orientation[4] * values_nd[19]
                H3 = self.axis_orientation[5] * values_nd[20]

        #---Pick the element----------------------------------------------------
            elem = self.fem_mesh.elements[member_id-1] #....................... Pick element at index

        #---Save linear and angular displacement of element---------------------
            elem.u1 = elem.u1 + self.relax_factor * (u1 - elem.u1)
            elem.u2 = elem.u2 + self.relax_factor * (u2 - elem.u2)
            elem.u3 = elem.u3 + self.relax_factor * (u3 - elem.u3)
            elem.t1 = elem.t1 + self.relax_factor * (t1 - elem.t1)
            elem.t2 = elem.t2 + self.relax_factor * (t2 - elem.t2)
            elem.t3 = elem.t3 + self.relax_factor * (t3 - elem.t3)

        #---Get stiffness and mass matrix at midpoint---------------------------
            stiff_mat_elem = (elem.stiffmat1 + elem.stiffmat2)/2 #............. Stiffness matrix at mid point
            mass_mat_elem  = (elem.massmat1  + elem.massmat1)/2 #.............. Stiffness matrix at mid point

        #---Get strain of the element-------------------------------------------
            strain_vec = np.dot(stiff_mat_elem, [F1, F2, F3, M1, M2, M3])

            elem.gamma11 = strain_vec[0] #..................................... beam axial stretching strain
            elem.gamma22 = strain_vec[1]/2 #................................... engineering transverse shear strains along x2
            elem.gamma13 = strain_vec[2]/2 #................................... engineering transverse shear strains along x3
            elem.k1      = strain_vec[3] #..................................... Twist
            elem.k2      = strain_vec[4] #..................................... curvature around x2
            elem.k3      = strain_vec[5] #..................................... curvature around x3

        #---Get velocities of element-------------------------------------------
            if self.analysis_type == "unsteady":
                inv_mass_mat_elem = np.linalg.inv(mass_mat_elem) #.............. Inverse mass matrix
                vel_vec = np.dot(inv_mass_mat_elem, [P1, P2, P3, H1, H2, H3]) #. Compute velocity vector
            #---Set values------------------------------------------------------
                elem.v1     = self.relax_factor * vel_vec[0]
                elem.v2     = self.relax_factor * vel_vec[1]
                elem.v3     = self.relax_factor * vel_vec[2]
                elem.theta1 = self.relax_factor * vel_vec[3]
                elem.theta2 = self.relax_factor * vel_vec[4]
                elem.theta3 = self.relax_factor * vel_vec[5]

    #---Return------------------------------------------------------------------
        return

    """
    ----------------------------------------------------------------------------
        Write input file
    ----------------------------------------------------------------------------
    """
    def __write_input_file(self):
        analysis_flag_type = {"steady":1, "unsteady":2} #...................... Flags of gebt analysis type

        file_path = self.file_path_input
        self.message.out("Write gebt input file for time step to: %s\n" %(file_path))

    #---Compute gain factor-----------------------------------------------------
        gain_factor=1.0
        if self.time_step < self.delay_iter:
            gain_factor = 0.0
        if self.ramp_step:
            if self.time_step >= self.delay_iter and self.time_step <= (self.delay_iter+self.ramp_step):
                gain_factor = 1.0 / self.ramp_step * (self.time_step - self.delay_iter)
        print("gain factor ", gain_factor)

        with open(file_path, "w") as f:
        #---Write analysis control parameters-----------------------------------
            analysis_flag = analysis_flag_type[self.analysis_type]
            niter = self.niter
            nstep = 100 #self.nstep

            f.write("%i " %(analysis_flag) )
            f.write("%i " %(niter)         )
            if analysis_flag < 1:
                f.write("%i " %(nstep))
            else:
                f.write("1 ")
            f.write("# analysis_flag niter nstep\n")
            f.write("\n"                   )

        #---Write motion parameters of the global frame of reference------------
            if analysis_flag > 0:
                f.write("0 0 0\n")
                f.write("0 0 0\n")
                f.write("0 0 0\n")
                f.write("0 0 0\n")
                f.write("\n"           )

        #---Next definition of geometry and boundary condition parameters-------
            nkp      = len(self.fem_mesh.nodes) #.............................. the total number of key points
            nmemb    = len(self.fem_mesh.elements)  #.......................... the total number of members
            ncond_pt = len(self.fem_mesh.nodes)+1  #........................... the total number of points having prescribed conditions (including boundary conditions) + 1 fixed node in space
            nmate    = 2*len(self.fem_mesh.elements)  #........................ the total number of cross-sections
            nframe   = len(self.fem_mesh.elements)  #.......................... the total number of frames
            ncond_mb = 0  #.................................................... the total number of members having prescribed, distributed loadings
            ndistr   = 0  #.................................................... the total number of functions used to approximate the distributed loads
            ntimefun = 0  #.................................................... the total number of time functions
            ncurv    = 0  #.................................................... total number of initial curvature/twist sets including k1 , k2 , k3

            f.write("%i " %(nkp)      )
            f.write("%i " %(nmemb)    )
            f.write("%i " %(ncond_pt-1) )
            f.write("%i " %(nmate)    )
            f.write("%i " %(nframe)   )
            f.write("%i " %(ncond_mb) )
            f.write("%i " %(ndistr)   )
            f.write("%i " %(ntimefun) )
            f.write("%i " %(ncurv)    )
            f.write("# nkp nmemb ncond pt nmate nframe ncond mb ndistr ntimefun ncurv" )
            f.write("\n")

        #---Print out list of key points----------------------------------------
            for k, kp in enumerate(self.fem_mesh.nodes):
                f.write("%i " %(kp.ID))
                f.write("%f " %(kp.x0))
                f.write("%f " %(kp.y0))
                f.write("%f " %(kp.z0))
                f.write("\n")
            f.write("\n")

        #---Write out element definition----------------------------------------
        #   (memb_no kp_1 kp_2 mate_no1 mate_no2 frame_no ndiv curv_no)
            for k, elem in enumerate(self.fem_mesh.elements):
                f.write("%i " %(elem.ID)) #.................................... memb_no
                f.write("%i " %(elem.node1.ID)) #.............................. kp_1
                f.write("%i " %(elem.node2.ID)) #.............................. kp_2
                f.write("%i " %(2*k+1)) #...................................... mate_no1
                f.write("%i " %(2*k+2)) #...................................... mate_no2
                f.write("%i " %(k+1)) #........................................ frame_no
                f.write("1 " ) #............................................... ndiv
                f.write("0 " ) #............................................... curv_no (unused, i.e. 0)
                f.write("\n")
            f.write("\n")

        #---Write point condition for fixed node--------------------------------
            fixed_nd = self.fem_mesh.fixed_nd #................................ Node, which is fixed in space
            f.write("%i\n" %(fixed_nd.ID))
            f.write("1 2 3 4 5 6\n")
            f.write("0 0 0 0 0 0\n")
            f.write("0 0 0 0 0 0\n")
            f.write("0 0 0 0 0 0\n")
            f.write("\n")

        #---Write point conditions for force and moments------------------------
            #f.write("#Forces and moments\n")
            for k, kp in enumerate(self.fem_mesh.nodes):
                if not kp == self.fem_mesh.fixed_nd:
                    f.write("%i\n" %(kp.ID))
                    f.write("7 8 9 10 11 12\n")
                    #f.write("%f %f %f %f %f %f\n" %(kp.fx, kp.fy, kp.fz, kp.mx, kp.my, kp.mz))
                    f.write("%e " %(gain_factor * kp.fx))
                    f.write("%e " %(gain_factor * kp.fy))
                    f.write("%e " %(gain_factor * kp.fz))
                    f.write("%e " %(gain_factor * kp.mx))
                    f.write("%e " %(gain_factor * kp.my))
                    f.write("%e " %(gain_factor * kp.mz))
                    f.write("\n"                        )

                    f.write("0 0 0 0 0 0\n") #................................. tf_1, tf_2, tf_3, tf_4, tf_5, tf_6
                    f.write("0 0 0 0 0 0\n") #................................. ff_1, ff_2, ff_3, ff_4, ff_5, ff_6
                    f.write("\n")
            #f.write("\n")

        #---Write out stiffness and mass matrix---------------------------------
            #f.write("#Stiffness and mass matrices\n")
            for k, elem in enumerate(self.fem_mesh.elements):
            #---First node------------------------------------------------------
                f.write("%i\n" %(2*elem.ID-1))
                self.__write_mat_to_file(f, elem.stiffmat1) #.................. Stiffness matrix for node 1
                if analysis_flag > 0:
                    self.__write_mat_to_file(f, elem.massmat1) #............... Mass matrix for node 1

            #---Second node-----------------------------------------------------
                f.write("%i\n" %(2*elem.ID))
                self.__write_mat_to_file(f, elem.stiffmat2) #.................. Stiffness matrix for node 2
                if analysis_flag > 0:
                    self.__write_mat_to_file(f, elem.massmat2) #............... Mass matrix for node 2

        #---Write out frames----------------------------------------------------
            #f.write("#Local frame of nodes\n")
            for k, elem in enumerate(self.fem_mesh.elements):
            #---Get directional frame for first node in this element------------
                frame_mat = elem.frame_matrix
            #...Write frame to file.............................................
                f.write("%i\n" %(elem.ID))
                self.__write_mat_to_file(f, frame_mat)

        #---Write simulation range----------------------------------------------
            if self.analysis_type == "unsteady":
                t1 = self.dt * self.time_step
                t2 = self.dt * (self.time_step+1)
                f.write("%e %e\n" %(t1, t2))

    #---Return------------------------------------------------------------------
        return

    """
    ----------------------------------------------------------------------------
        Write ini file (only unsteady mode)
    ----------------------------------------------------------------------------
    """
    def __write_ini_file(self):
        path_file_ini = self.file_path_input + ".ini"
        with open(path_file_ini, "w") as f:
            for k in range(self.ini_conditions.shape[0]):
                values = self.ini_conditions[k,:]
                f.write("%e " %(self.ini_conditions[k,0]))
                f.write("%e " %(self.ini_conditions[k,1]))
                f.write("%e " %(self.ini_conditions[k,2]))
                f.write("%e " %(self.ini_conditions[k,3]))
                f.write("%e " %(self.ini_conditions[k,4]))
                f.write("%e " %(self.ini_conditions[k,5]))
                f.write("\n")
    #---Return------------------------------------------------------------------
        return

    """
    ----------------------------------------------------------------------------
        Write matrix to file
    ----------------------------------------------------------------------------
    """
    def __write_mat_to_file(self, f, mat):
        Ni = mat.shape[0]
        Nj = mat.shape[1]
        for i in range(Ni):
            for j in range(Nj):
                f.write("%0.16e " %(mat[i,j]))
            f.write("\n")
        f.write("\n")
        return

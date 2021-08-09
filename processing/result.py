import os
import numpy as np

class variable:
    """
    attributes:
        name  - name
        dof   - degrees of freedom
        np    - number of points
        idx_0 - starting index from permuation table

    """
    def __parse_varname(self, line):
        if "[" in line:
            name = line.split("[")[0]
        else:
            name = line.split(":")[0].strip()
        return name

    def __parse_vardof(self, line):
        if "[" in line:
            dof = line.split("]")[1].split(":")[1].strip().split()[-1].strip()
        else:
            dof = line.split(":")[1].strip().split()[-1].strip()
        return int(dof)

    def __parse_np(self, line):
        """ find the number of points per variable this is necessary becuase the
        free-surface is solved for on it's own body, a lower dimensional body with
        less points than the rest of the mesh
        """
        if "[" in line:
            var_np = line.split("]")[1].split(":")[1].strip().split()[-1].strip()
        else:
            var_np = line.split(":")[1].strip().split()[0].strip()
        return int(var_np)

    def __init__(self, line):
        self.name  = self.__parse_varname(line)
        self.dof   = self.__parse_vardof(line)
        self.np    = self.__parse_np(line)

    def parse_idx0(self, line):
        # TO DO: parse idx0 from permuation table
        pass

class result_file:
    def __parse_nodes(self):
        """Extract coordinate and node_id's from mesh.nodes file
        """
        with open(self.nodes_fp, "r") as file:
            nodes = file.readlines()
            n_nodes = len(nodes)

            coords = np.empty(
                n_nodes,
                dtype=[
                    (
                        "node_id",
                        np.float64,
                    ),  # this shouldn't be float64 but needs to be for np.hstack
                    ("x", np.float64),
                    ("y", np.float64),
                    ("z", np.float64),
                ],
            )

            for i in range(n_nodes):
                node_id, _, x, y, z = nodes[i].split()
                coords[i] = node_id, x, y, z

            # Not sure how important this is?
            del nodes

        # Set coords array as an attribute
        self.coords = coords.copy()

    def __parse_mesh_params(self, simul_type="Steady"):
        """Find mesh parameters (ie. Nx, Ny, Nt) from mesh.nodes and .result
        """

        mesh_params = {
            "Nx": len(np.unique(self.coords["x"])),
            "Ny": len(np.unique(self.coords["y"])),
            "Nnodes": len(np.unique(self.coords["node_id"])),
            "Simulation Type": simul_type.lower(),
        }

        # This might be a faster way to do this:
        #       https://stackoverflow.com/a/35857833/10221482
        # or this, since we are only couting one word:
        #       https://stackoverflow.com/a/38401151/10221482

        # Sweep .result file and find Nt by counting occurance of the string "Time".
        with open(self.result_fp, "r") as file:
            # count the number of ocurrance of "Time" (i.e.) number of timesteps
            nt = sum(map(lambda x: "Time:" in x, file.readlines()))

        if mesh_params["Simulation Type"] == "steady":
            mesh_params["S.S. itters"] = nt
        else:
            raise NotImplemented("time dependent data not yet supported")
        self.mesh_params = mesh_params

    def __read_result(self):
        # open the result file
        data_file = open(self.result_fp, "r")

        nt = self.mesh_params["S.S. itters"]
        nn = self.mesh_params["Nnodes"]

        ########################################################################
        # Parse the header. For more info see (pg. 127):
        # http://www.nic.funet.fi/pub/sci/physics/elmer/doc/ElmerSolverManual.pdf
        ########################################################################
        while True:
            line = data_file.readline()
            if "Degrees of freedom:" in line:
                var_start = data_file.tell()
            if "Total DOFs:" in line:
                var_end = data_file.tell()
            if "Number Of Nodes:" in line:
                header_end = data_file.tell()
                break

        ########################################################################
        # Parse vairable names and dofs from header. For info see above ^^^
        ########################################################################
        data_file.seek(var_start)
        variables = []
        while True:
            line = data_file.readline()
            if "Total DOFs:" not in line:
                # creat instance of variable class
                var = variable(line)
                if int(var.dof) < 2:
                    variables.append(var)
            else:
                break

        ########################################################################
        # Allocate arrays to store the data
        ########################################################################
        dstart = np.zeros([len(variables), nt], dtype=np.int)
        count = np.zeros(len(variables), dtype=np.int)
        data = np.empty([nn, nt], dtype=[(var.name, np.float64) for var in variables])

        # # loop through header one more time to find the number of points per variable
        # # this is necessary becuase the free-surface is solved for on it's own
        # # body, a lower dimensional body with less points than the rest of the mesh
        # data_file.seek(var_start)
        # while True:
        #     line = data_file.readline()
        #     if "Total DOFs:" not in line:
        #         for i, variable in enumerate(variables):
        #             if variable in line:
        #                 var_np[i] = self.__parse_np(line)
        #     else:
        #         break
        ########################################################################
        # Read data for each variable for each timestep (or itteration)
        ########################################################################
        # Only start scanning after the header
        data_file.seek(header_end)

        # Find the postion of variable of interest
        while True:
            line = data_file.readline()

            # Once at the end of the file, break
            if not line:
                break

            for k, var in enumerate(variables):
                if var.name + "\n" == line:
                    i = count[k]
                    dstart[k, i] = data_file.tell()
                    count[k] += 1
                if count[k] > nt:
                    raise IndexError(
                        "There were more matches than timesteps. Somethings wrong"
                    )

        # Go to the postion in the .result file
        for k, var in enumerate(variables):
            for i in range(nt):
                data_file.seek(header_end)
                data_file.seek(dstart[k, i])

                line = data_file.readline()
                if line[6:] != "use previous\n":
                    for j in range(var.np):
                        data_file.readline()

                for j in range(var.np):
                    data[var.name][j, i] = float(data_file.readline())

        # All done!
        data_file.close()

        # Only returing the final itteration since
        self.data = data[:, -1]

    def __init__(self, result_fp, nodes_fp, simul_type="Steady"):

        self.result_fp = result_fp
        self.nodes_fp = nodes_fp

        # parse the nodes file
        self.__parse_nodes()
        # parse the mesh paramters
        self.__parse_mesh_params(simul_type="Steady")
        # parse the data from the results file
        self.__read_result()

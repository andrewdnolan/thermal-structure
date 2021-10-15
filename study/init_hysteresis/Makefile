# folder to write the mesh too
MESH_DIR = mesh
# name of the .grd file to make the mesh form

# Mesh files to be made
mesh_extensions := boundary elements header nodes
# For each "word" in foo, output that same word with an exclamation after
mesh_files = $(foreach wrd,$(mesh_extensions),$(MESH_DIR)/mesh.$(wrd))


all: $(mesh_files)

# make the mesh dir if does not exist
$(mesh_files) : $(MESH_DIR)
	cd  ../.. && ./scripts/make_mesh.sh 200 ./study/init_hysteresis/$(MESH_DIR)

# make the mesh dir if does not exist
$(MESH_DIR) :
	test ! -d $(MESH_DIR) && mkdir $(MESH_DIR)

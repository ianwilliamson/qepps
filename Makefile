EXE=qepps
SRC_DIR = src

.PHONY: qepps

$(EXE):
	$(MAKE) -C $(SRC_DIR) $(EXE)

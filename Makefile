SRC_DIR=src
LUA_DIR=$(SRC_DIR)/lua

all: lua
	$(MAKE) -C $(SRC_DIR) all
	
lua:
	$(MAKE) -C $(LUA_DIR)


.PHONY: clean clean-lua

clean:
	rm qepps
	$(MAKE) -C $(SRC_DIR) clean

clean-lua:
	$(MAKE) -C $(LUA_DIR) clean

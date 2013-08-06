int create_mapping();

int close_mapping();

int write_to_shared_memory(PARA_DATA *para, REAL **var);

int read_from_shared_memory(PARA_DATA *para, REAL **var);

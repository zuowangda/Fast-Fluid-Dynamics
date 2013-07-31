int create_shared_memory(char *ffd_memory_name, char *other_memory_name);

int create_mapping();

int close_mapping();

int write_to_shared_memory(ffdSharedData *ffdData);

int read_from_shared_memory(otherSharedData *otherData);


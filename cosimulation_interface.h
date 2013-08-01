int create_mapping(char *ffdDatNam, char *modDatNam);

int close_mapping();

int write_to_shared_memory(ffdSharedData *ffdData);

int read_from_shared_memory(otherSharedData *otherData);


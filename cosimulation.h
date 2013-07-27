

int read_cosimulation_data(PARA_DATA *para, REAL **var);

void receive_data();

BOOL CALLBACK receive_data_dialog(HWND hwndDlg, UINT message, 
                                  WPARAM wParam, LPARAM lParam);

BOOL CALLBACK receive_command(HWND hwndDlg, UINT message, WPARAM wParam, 
                              LPARAM lParam);

void send_command(SentCommand command_sent);

int send_data(SentData data_sent);

int write_cosimulation_data(PARA_DATA *para, REAL **var);

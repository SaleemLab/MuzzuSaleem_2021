function users = getUserList
% This just loads the user list.
% AS 09/4/2019 
DIRS = SetDefaults;
users  = readtable([DIRS.DatabaseUsers 'UserList.csv']);
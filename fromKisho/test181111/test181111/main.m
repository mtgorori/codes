close all;
clear;

logger('-------------------------------------\n');
logger('This is test program of logger\n');
logger('date: %s\n', datestr(datetime));
logger('-------------------------------------\n');

numint = 23;
numdouble = 174.52;
line = 'John Lennon';


logger('Nice to meet you\n');
logger('My name is %s.\n', line);
logger('I am %d years old and %f cm tall\n', numint, numdouble);
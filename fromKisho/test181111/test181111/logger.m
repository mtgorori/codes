function y = logger(formatSpec, varargin)
% ログ用の関数: 標準出力およびログファイルへの出力を行う
% 基本的には標準出力のfprintfと同じように使えばよい
% ログには呼び出し元の関数名と呼び出した行数が入る
% この関数は呼び出し毎にファイルを開いている

% formatSpec: 書式文字列
% varargin: 書式に入れ込む変数（複数）

st = dbstack(); % 関数スタックに関する構造体
callfname = st(2).name; % 呼び出し元の関数の名前
callline = num2str(st(2).line); % 呼び出し元関数での呼び出し行

% 関数の名前と呼び出し行は文字列で与えられる
% これらをログ上で並べたとき，きれいに見えるように間隔調整をする
len1 = strlength(callfname); 
len2 = strlength(callline);
maxlen1 = 17;
maxlen2 = 4;
bl1 = blanks(abs(maxlen1 - len1)); % 関数名の後の空白の数
bl2 = blanks(abs(maxlen2-len2)); % 行数名の後の空白の数

% 空白による各文字列の間隔調整
formatSpec2 = ['$ [' callfname '] line:' callline bl1 bl2 '| ' formatSpec]; 

%-----------------
% 標準出力に出す
%-----------------

fprintf(formatSpec2, varargin{1:nargin-1}); 

%------------------
% ファイルへの出力
%------------------

fid = fopen('./logger.txt', 'a'); % ファイルを開く

fprintf(fid, formatSpec2, varargin{1:nargin-1}); % ファイルに書き込み
fclose(fid); % ファイルを閉じる

y = 0;
end
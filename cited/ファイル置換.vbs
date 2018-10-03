Option Explicit

Const TARGET_FOLDER = "H:\backup\codes\immature"     '※フォルダを指定
Const REPLACE_FROM = "\\Azlab-fs01\東研究室\個人work\竹内(ひ)" '※
Const REPLACE_TO   = "H:\backup\" '※

Const ForReading = 1 '読み込み
Const ForWriting = 2 '書きこみ（上書きモード）
Const ForAppending = 8 '書きこみ（追記モード）

Dim objFSO, objFolder, objFile, objSubFolder, objTXT
Set objFSO = WScript.CreateObject("Scripting.FileSystemObject")
Set objFolder = objFSO.GetFolder(TARGET_FOLDER)

For Each objFile In objFolder.Files
    Dim strFilePath, infile, outfile, strData
    strFilePath = objFSO.BuildPath(TARGET_FOLDER, objFile.Name)
    Set infile = objFSO.OpenTextFile(strFilePath,ForReading)
    strData = infile.ReadAll
    infile.Close
    Set infile = Nothing
    Set outfile = objFSO.OpenTextFile(strFilePath,ForWriting)   '(上書き)
    outfile.Write Replace(strData,REPLACE_FROM,REPLACE_TO)
    outfile.Close
    Set outfile = Nothing
Next

Set objFolder = Nothing
Set objFSO = Nothing

MsgBox "終了", vbInformation

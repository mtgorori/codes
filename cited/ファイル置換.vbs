Option Explicit

Const TARGET_FOLDER = "H:\backup\codes\immature"     '���t�H���_���w��
Const REPLACE_FROM = "\\Azlab-fs01\��������\�lwork\�|��(��)" '��
Const REPLACE_TO   = "H:\backup\" '��

Const ForReading = 1 '�ǂݍ���
Const ForWriting = 2 '�������݁i�㏑�����[�h�j
Const ForAppending = 8 '�������݁i�ǋL���[�h�j

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
    Set outfile = objFSO.OpenTextFile(strFilePath,ForWriting)   '(�㏑��)
    outfile.Write Replace(strData,REPLACE_FROM,REPLACE_TO)
    outfile.Close
    Set outfile = Nothing
Next

Set objFolder = Nothing
Set objFSO = Nothing

MsgBox "�I��", vbInformation

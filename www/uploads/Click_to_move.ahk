CheckItemLevel(){
BlockInput On
Send {LButton}
Send {Enter}
Sleep 2
Send /itemlevel
Send {Enter}
Sleep, 75
Send {LButton}
BlockInput Off
return
}
F1::CheckItemLevel() ; Assign F10 to Check item level
;----Remaining----	;Shows how many monsters are remaining in the instance.
F2::
Send {Enter}/remaining{Enter}
return

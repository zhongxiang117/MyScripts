# check all `ps1` files:
#$PROFILE | Select-Object *

# list all environments
#ls env:* | Sort-Object Name


# to create a file
#if (!(Test-Path -Path $PROFILE)) {New-Item -ItemType File -Path $PROFILE -Force}


# to reload settings
# . \path\to\profile.ps1


# color theme -> `Windows-Terminal` -> Settings -> JSON
#
#{
#    "name": "Ubuntu",
#    "black": "#2e3436",
#    "red": "#cc0000",
#    "green": "#4e9a06",
#    "yellow": "#c4a000",
#    "blue": "#3465a4",
#    "purple": "#75507b",
#    "cyan": "#06989a",
#    "white": "#d3d7cf",
#    "brightBlack": "#555753",
#    "brightRed": "#ef2929",
#    "brightGreen": "#8ae234",
#    "brightYellow": "#fce94f",
#    "brightBlue": "#729fcf",
#    "brightPurple": "#ad7fa8",
#    "brightCyan": "#34e2e2",
#    "brightWhite": "#eeeeec",
#    "background": "#300a24",
#    "foreground": "#eeeeec"
#}


# set terminal `bash`-like
Set-PSReadLineOption -EditMode 'Emacs'
Set-PSReadLineKeyHandler -Key Tab -Function MenuComplete
Set-PSReadlineOption -BellStyle None

# for miniconda
# !! Contents within this block are managed by 'conda init' !!
If (Test-Path "C:\Users\xiang\miniconda3\Scripts\conda.exe") {
    (& "C:\Users\xiang\miniconda3\Scripts\conda.exe" "shell.powershell" "hook") | Out-String | ?{$_} | Invoke-Expression
}

# it has to be put after miniconda
function prompt {
    $host.ui.RawUI.WindowTitle = "$pwd"
    $curdir = ($pwd).Path.Split('\')
    if ($curdir[-1] -eq '') {
        $curdir = $curdir[0] + "\\"
    } else {
        $curdir = $curdir[-1] + ":\\"
    }
    $cudaprefix = $env:CONDA_PROMPT_MODIFIER
    if ($cudaprefix -ne '') {
        Write-Host "${cudaprefix}level-$env:CONDA_SHLVL @ " -NoNewline
    }
    Write-Host "$curdir" -ForegroundColor Green -NoNewline
    "> "
}
# (Get-Command Prompt).ScriptBlock




# bash

```sh
setxkbmap -option caps:escape

alias gdb="gdb --silent"
```

# vim
```vim
" show line numbers (nu)
set number

" display relative line numbers instead of absolute (rnu)
set relativenumber

" the width of indent for operations like << and >> (sw)
set shiftwidth=4

" display size of the tab character (ts)
set tabstop=4

" use the indentation of the current line when creating a new one (ai)
set autoindent

" constantly show the file name (ls)
set laststatus=2

"" not used

" set vim colorscheme (colo)
colorscheme delek

" autoindent when typing { and }
set smartindent

" number of spaces that a tab counts when typing <TAB> (sts)
set softtabstop=4

" inserted tabs turn into spaces (et)
set expandtab

" show information about the current position (activated by default, but just
" in case)
set ruler

" show tabs and trailing spaces
set list
```
_fiesta()
{
    _init_completion || return

    _filedir '@(lua)'
}

complete -F _fiesta fiesta

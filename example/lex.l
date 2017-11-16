EMPTY \s\s*
BREAK ;
RESERVED var|function|return|if|else|for|while
OPERATOR [\(\)\{\}\[\]\+\-\*/=(\+=)(\-=)(\*=)(/=)><(>=)(<=)(!=)(\+\+)(\-\-)!&\|(&&)(\|\|)(&=)(\|=)]
ID [a-zA-Z_][a-zA-Z_0-9]*
NUMBER [0-9][0-9]*
REAL [0-9][0-9]*.[0-9][0-9]*
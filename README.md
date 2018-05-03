# Simple-Lex

A very simple version of Lex that generates lexical analyzers.

It reads an input stream specifying the lexical analyzer and outputs source code implementing the lexer in the C++ programming language.

For more information about the real Lex, see [Lex (software)](https://en.wikipedia.org/wiki/Lex_(software)).

## Goal

I program it to help myself understand the principle of the *Lex* better.

## Simplification

The *.l* file can only include token names and their regular expressions that indicate the rules.

## Example

The *lex.l* is as follow.

```
EMPTY \s\s*
BREAK ;
RESERVED var|function|return|if|else|for|while
OPERATOR [\(\)\{\}\[\]\+\-\*/=(\+=)(\-=)(\*=)(/=)><(>=)(<=)(!=)(\+\+)(\-\-)!&\|(&&)(\|\|)(&=)(\|=)]
ID [a-zA-Z_][a-zA-Z_0-9]*
NUMBER [0-9][0-9]*
REAL [0-9][0-9]*.[0-9][0-9]*
```

Execute `lex lex.l`.

And a source file named *lex.yy.cpp* will be generated.

Compile it: `g++ lex.yy.cpp -o a`.

The code file named *src.txt* that is required to be analyzed:

```javascript
var sum = 0;
var real = 1.1;
function count(arg) {
    for (var i = 0; i < arg; ++i) {
        sum += 1;
        real -= 2.2;
    }
}
count(10);
```

Execute `a src.txt`.

The result is such a token-stream printed in the console, whose pattern is `<token_name, token_word>` such as `<RESERVED, "var">`.

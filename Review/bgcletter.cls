%%
%% This is file `bgcletter.cls',
%% modified from the mpslet.cls created by Patrick W. Daly
%% 
%% 
\NeedsTeXFormat{LaTeX2e}[1994/12/01]
\ProvidesClass{bgcletter}
          [2011/01/19 1.3a (PWD)]
 %-------------------------------------------------------------------
 % NOTICE:
 % This file may only be used for producing official letters
 %  of the Max-Planck-Institut f\"ur Biogeochemie
 %
%% Copyright (C) 2009 Max-Planck-Institut f\"ur Biogeochemie
 % Hans-Kn"\oll-Str. 10
 % 07745 Jena
 % Germany
 %
 %  Written by Carlos A. Sierra, 
 %  based on earlier versions by Patrick W. Daly,
 %  and for MPI Aeronomie (MPAe), predecessor to MPS
 %
 %-----------------------------------------------------------
 % A LaTeX Class for official MPI-BGC letters
 % Call as  \documentclass[<options>]{bgcletter}
 % where <options> can be
 %     10pt  11pt   12pt  for basic font size (default=10pt)
 %     colour | color     for colour output (def)
 %     bw | sw            for black&white output
 %     nologo             for no logos, intended to print on preprinted paper
 %     logo               to output MPG header and footer on first page (def)
 %     doc                adds name of source file to head info
 %     nodoc              suppresses adding name of source file (def)
 %     german             loads german.sty and translates some texts to German
 %     english            selects (UK) English, ie English texts and UK date
 %
 % Enter letter information with
 %   \name{...} \signature{...} (if \signature omitted, then name used)
 %   \telephone{..} \internet{..} \email{..} (Optional)
 %   \fax{..}  \workstation{..} (last defaults to mps)
 %   \homephone{..} \homefax{..} \mobile{..} (with Vorwahl without 0)
 %   \myref{..} \yref{..} \ymail{..} \subject{..} \date{..} (all optional)
 %
 %-----------------------------------------------------------


\newif\if@german \@germanfalse
\DeclareOption{german}{\@germantrue}
\DeclareOption{ngerman}{\@germantrue}
\DeclareOption{english}{\@germanfalse}
\newif\if@kolor \@kolortrue
\DeclareOption{colour}{\@kolortrue}
\DeclareOption{color}{\@kolortrue}
\DeclareOption{bw}{\@kolorfalse}
\DeclareOption{sw}{\@kolorfalse}
\newif\if@logo \@logotrue
\DeclareOption{nologo}{\@logofalse\@kolorfalse}
\DeclareOption{logo}{\@logotrue}
\newif\if@docname \@docnamefalse
\DeclareOption{nodoc}{\@docnamefalse}
\DeclareOption{doc}{\@docnametrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}
\ProcessOptions
\LoadClass[a4paper]{letter}
\RequirePackage{ngerman}
\if@german\else\selectlanguage{english}\fi
\if@kolor
  \RequirePackage{color}
  \definecolor{mpg}{rgb}{.187,.613,.594}
  \newcommand*\mpgcol{\leavevmode\color{mpg}}
  \newcommand*\bgclogo{bgc_color}
  \newcommand*\mpglogo{minerva}
\else
  \let\mpgcol=\relax
  \newcommand*\bgclogo{mps_grau}
  \newcommand*\mpglogo{minerva_logo_bw}
\fi
\RequirePackage{graphicx}
\newcommand{\titlea}{\fontsize{19}{24}\usefont{\encodingdefault}{pnc}{m}{n}%
  \mpgcol}
\newcommand{\titleb}{\fontsize{10.6}{13}\usefont{\encodingdefault}{pnc}{m}{n}%
  \mpgcol}
\newcommand{\addfont}{\fontsize{7}{8}\usefont{\encodingdefault}{phv}{m}{n}%
  \mpgcol}
\newcommand{\senderfont}{\fontsize{8}{9}\usefont{\encodingdefault}{\rmdefault}{m}{n}}
\def\gernames{\protect\typeout{German Names}}
\def\engnames{\protect\typeout{English Names}}
\newcommand{\defGerEng}[3]{%
  \let\Def\relax
  \let\protect\@unexpandable@protect
  \edef\gernames{\gernames\Def{#1}{#2}}%
  \edef\engnames{\engnames\Def{#1}{#3}}
  \let\Def\@namedef
  \let\protect\relax}
\newcommand*{\MP@names}{\ifnum\language=\l@ngerman\gernames
     \else\engnames\fi}
\AtBeginDocument{\providecommand{\l@ngerman}{1}}
\renewcommand*{\name}[2][\@empty]
    {\ifx\@empty#1\defGerEng{fromname}{#2}{#2}\else
        \defGerEng{fromname}{#2}{#1}\fi}
\renewcommand*{\signature}[2][\@empty]
    {\ifx\@empty#1\defGerEng{fromsig}{#2}{#2}\else
        \defGerEng{fromsig}{#2}{#1}\fi}
\defGerEng{yrefname}    {Ihr Zeichen:}{Your Ref.:}
\defGerEng{myrefname}    {Unser Zeichen:}{Our Ref.:}
\defGerEng{ymailname}   {Ihr Schreiben vom}{Your letter of}
\defGerEng{subjectname} {\textsl{Betr.}\,}{\textsl{Subject}}
\defGerEng{docname}     {Dok:}{Doc:}
\defGerEng{enclname}{Anlagen}{encl}
\defGerEng{ccname}{Kopien an}{cc}
\defGerEng{phone}{(0)3641-57-}{+49-(0)3641-57-}
\defGerEng{intcode}{0}{+49-}
\newcommand*\mydomain{.mpg.de}
\newcommand*\domain{@bgc-jena.mpg.de}
\newcommand*\webpage{www.bgc-jena.mpg.de}
\newcommand*\precode{\,+\,49-(0)\,3641~/\,57\,--\,}
\newcommand*\firmgerm{Max--Planck--Institut f\"ur Biogeochemie}
\newcommand*\firmeng{Max Planck Institute for Biogeochemistry}
\newcommand*\firmadd{MPI f\"ur Biogeochemie $\cdot$ Postfach 10 \, 01 \, 64
      $\cdot$ 07745 Jena, Germany}
\InputIfFileExists{bgcletter.def}{}{}
\InputIfFileExists{bgcletter.cfg}{\typeout{%
*******************************************^^J%
* Read in local configuration bgcletter.cfg *^^J%
*******************************************}}{}
\setlength{\oddsidemargin}  {-5mm} % +1in = 20mm
\setlength{\evensidemargin} {-20mm}% +1in = 5mm
\setlength{\marginparwidth} {90pt}
\setlength{\marginparsep}   {11pt}
\setlength{\marginparpush}  {5pt}
\setlength{\topmargin}      {-20mm}
\setlength{\headheight}     {16mm}

\setlength{\headsep}        {25pt}
\setlength{\footskip}       {30pt}
\setlength{\textheight}     {\paperheight}
 \addtolength{\textheight}{-1.75in}
 \addtolength{\textheight}{-\topmargin}
 \addtolength{\textheight}{-\headheight}
 \addtolength{\textheight}{-\headsep}
 \addtolength{\textheight}{-\footskip}

\setlength{\textwidth}{137mm}
\newlength{\fullwidth}
\setlength{\fullwidth}{\paperwidth}
 \addtolength{\fullwidth}{-2in}
 \addtolength{\fullwidth}{-\oddsidemargin}
 \addtolength{\fullwidth}{-\evensidemargin}
\setlength{\columnsep}      {10pt}
\setlength{\columnseprule}  {0pt}
\raggedbottom
\setlength{\longindentation}{.4\textwidth}
\setlength{\indentedwidth}{\textwidth}
\addtolength{\indentedwidth}{-\longindentation}
\newlength{\leftfield}
\setlength{\leftfield}{\textwidth}
\addtolength{\leftfield}{5\p@}
\newlength{\rightfield}
\setlength{\rightfield}{\fullwidth}
\addtolength{\rightfield}{-\leftfield}
\let\MPletter=\letter
\def\MP@letter#1\@?@?{\MPletter{#1}}
\renewcommand{\letter}[2][]{\gdef\alttoname{#1}%
      \expandafter\MP@letter#2\@?@?}
\renewcommand{\opening}[1]{\MP@names\thispagestyle{firstpage}
 \vspace*{1.3cm}%
 \parbox[t][2.5cm][t]{\leftfield}{%
    \toname\\\toaddress}%
 \raisebox{-3mm}{\parbox[t][4cm][t]{\rightfield}{\senderfont
    \ifx\@empty\fromname\else\textbf{\fromname}\fi
    \ifx\@empty\telephonenum\else\\
       Tel.: \phone\telephonenum \fi
    \ifx\@empty\@homephone\else\\
       Pri.: \intcode\@homephone \fi
     \ifx\@empty\@fax\else\\
       Fax: \phone\@fax \fi
     \ifx\@empty\@homefax\else\\
       Fax: \intcode\@homefax \fi
     \ifx\@empty\@mobilenum\else\\
       Mob: \intcode\@mobilenum \fi
     \ifx\@empty\@internet\else\\
           \@internet @\@workstation\mydomain\fi
     \vspace*{1ex}
     \ifx\@empty\@myref\else\\
        \myrefname\ \@myref\fi
     \ifx\@empty\@yref\else\\
        \yrefname\ \@yref\fi
     \ifx\@empty\@ymail\else\\
        \ymailname\ \@ymail\fi
     \if@docname\\\docname\ \jobname.tex\fi
    }}
 \hspace*{\leftfield}%
 \makebox[\rightfield][l]{\@date}\par
         \ifx\@empty\@subject \else\@subject\par\fi
 \vspace{2\parskip}
 #1\par\nobreak}
\newcommand{\rclosing}[1]{\par\nobreak\vspace{\parskip}
 \stopbreaks
 \noindent
 \hspace*{\longindentation}%
 \parbox{\indentedwidth}{\centering
 \ignorespaces #1\\[6\medskipamount]
 \ifx\@empty\fromsig
   \ifx\@empty\fromname\else\fromname\fi
 \else \fromsig \fi\strut}
 \par}
\newcommand{\lclosing}[1]{\par\nobreak\vspace{\parskip}
 \stopbreaks
 \noindent
 \parbox{\indentedwidth}{\raggedright
 \ignorespaces #1\\[6\medskipamount]
 \ifx\@empty\fromsig
   \ifx\@empty\fromname\else\fromname\fi
 \else \fromsig \fi\strut}
 \par}

\let\closing\lclosing
\newcommand{\nosigclosing}[1]{\par\nobreak\vspace{\parskip}
 \stopbreaks
 \noindent
 \hspace*{\longindentation}%
 \parbox{\indentedwidth}{\centering \ignorespaces #1}
 \par}
\renewcommand{\ps@firstpage}{%
  \def\@oddhead{\if@logo\parbox[b][\z@][t]{\fullwidth}{%
    {\titlea \firmgerm}\\
    {\titleb \firmeng}\\[1.0cm]
    {\addfont \firmadd}\hspace{100mm}%
      \raisebox{-1mm}[\z@][\z@]{\makebox[\z@][r]
      {\includegraphics[width=5cm,clip]{\bgclogo}}\hfill}}\fi}
   \def\@oddfoot{\if@logo\firmfoot\fi}%
   \def\@oddfoot{\if@logo\firmfoot\fi}%
\def\@evenhead{}\def\@evenfoot{}}
\newcommand*{\firmfoot}{\makebox[\textwidth]{\addfont
  \begin{tabular}[t]{@{}l}
    Max-Planck-Institut f\"ur\\
    Biogeochemie\\
    Hans-Kn\"{o}ll-Stra{\ss}e 10\\
    07745 Jena\\
    Germany
  \end{tabular}
  \hfill
  \begin{tabular}[t]{@{}l@{}}
   Tel.:  \precode 6110\\
   Fax.:  \precode 6110\\
   http:// www.bgc-jena.mpg.de
  \end{tabular}
  \hfill
  \begin{tabular}[t]{@{}r@{}}
  Direktorium\\
  Markus Reichstein  (Managing Dir.)\\
  Martin Heimann\\
  Susan Trumbore\\
  ID-Nr. DE 129517720
  \end{tabular}
}%
  \hspace{20mm}\raisebox{-13mm}[\z@]{\makebox[\z@][l]{%
    \includegraphics[width=2cm]{\mpglogo}}}\hfill
  }
\renewcommand{\ps@headings}{%
  \def\@oddhead{\parbox{0.95\fullwidth}{{\titleb\firmgerm}\\[10pt]
     \slshape\headtoname\
       \ifx\alttoname\@empty\toname\else\alttoname\fi
       \hfill\@date\hfill\pagename\ \thepage\\
     }}
  \def\@oddfoot{}}
\if@logo
  \ps@headings
\else
  \ps@plain
\fi
\def\yref#1{\def\@yref{#1}}            \def\@yref{}
\def\ymail#1{\def\@ymail{#1}}         \def\@ymail{}
\def\myref#1{\def\@myref{#1}}         \def\@myref{}
\def\subject#1{\def\@subject{\textbf{#1}}}   \def\@subject{}
\def\internet#1{\def\@internet{\lowercase{#1}}}   \def\@internet{}
\let\email=\internet
\def\workstation#1{\def\@workstation{#1}}    \def\@workstation{bgc-jena}
\def\fax#1{\def\@fax{#1}}           \def\@fax{}

\def\homephone#1{\def\@homephone{#1}} \def\@homephone{}
\def\homefax#1{\def\@homefax{#1}} \def\@homefax{}
\def\mobile#1{\def\@mobilenum{#1}} \def\@mobilenum{}
%% End of file bgcletter.cls.
\endinput
%%
%% End of file `bgcletter.cls'.

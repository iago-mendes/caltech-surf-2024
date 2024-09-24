\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{document}[Document]
\LoadClass[12pt, twoside]{article}

% LIBRARIES
	\RequirePackage[
		inner=1.25in,
		outer=1in,
		top=1in,
		bottom=1in,
		headsep=0.5in,
		footskip=0.5in
	]{geometry}
	\RequirePackage[utf8]{inputenc}
	\RequirePackage{graphicx}
	\RequirePackage[labelfont=bf, width=.8\textwidth]{caption}
	\RequirePackage{float}
	\RequirePackage{amsmath, amsfonts, amssymb, esint, amsthm, dsfont}
	\RequirePackage{titlesec}
	\RequirePackage{color}
	\RequirePackage{hyperref}
	\RequirePackage{pgf,tikz,pgfplots, tikz-3dplot}
		\usetikzlibrary{arrows}
		\pgfplotsset{compat=1.15}
	\RequirePackage{mathrsfs}
	\RequirePackage{ifthen}
	\RequirePackage{enumitem}
	\RequirePackage{multicol}
	\RequirePackage{subcaption}
	\RequirePackage{fancyhdr}
	\RequirePackage[sorting=none]{biblatex}
	\RequirePackage[bottom]{footmisc}
	\RequirePackage{etoolbox}
	\RequirePackage{orcidlink}
	\RequirePackage{indentfirst}
	\RequirePackage{titling}


% SETUP
  % allow large figures in pages
	\renewcommand{\floatpagefraction}{.7}

	\definecolor{red}{rgb}{1,0,0}
	\definecolor{green}{rgb}{0,1,0}
	\definecolor{blue}{rgb}{0,0,1}
	\definecolor{link}{HTML}{387494}

	\hypersetup {
		colorlinks=true,
		linkcolor=black,
		filecolor=link,      
		urlcolor=link,
		citecolor=link,
	}

	\newcommand{\red}[1]{\textcolor{red}{#1}}
	\newcommand{\green}[1]{\textcolor{green}{#1}}
	\newcommand{\blue}[1]{\textcolor{blue}{#1}}
	\newcommand{\gray}[1]{\textcolor{gray}{#1}}


% HEADERS & FOOTERS
	\fancypagestyle{default} {
		\fancyhead{} % clear all header fields
		\renewcommand{\headrule}{\hrulefill}
		\fancyhead[OR]{\it \thepage}
		\fancyhead[EL]{\it \thepage}

		\fancyfoot{} % clear all footer fields
	}
	\pagestyle{default}
	\fancypagestyle{title} {
		\fancyhead{} % clear all header fields
		\renewcommand{\headrule}{}

		\fancyfoot{} % clear all footer fields
		\fancyfoot[C]{\it \thepage}
	}

% SECTIONS
	\titleformat
		{\section} % command
		{\bf\large} % shape
		{\thesection. } % label
		{0cm} % separation
		{} % before code
		[] % after code

	\titleformat
		{\subsection} % command
		{\it\large} % shape
		{\thesubsection. } % label
		{0cm} % separation
		{} % before code
		[] % after code


% COMMANDS

	\newcommand{\header}[3][noSubtitle] {
		\begin{center}
			\vspace*{0.5cm}

			% Title
			{\Large\bf #3}

			\vspace*{0.5cm}
			
			% Subtitle
			\def\argument{#1} \def\default{noSubtitle}
			\ifx\argument\default
				{}
			\else
				{\large #1}
			\fi

			\vspace*{0.1cm}

			\large
			
			% Author
			#2
			
			\vspace*{0.1cm}

			% Date
			\today

			\vspace*{1cm}
		\end{center}
	}

	\newcommand{\instructions}[1]{
		\gray{#1} \\
	}

	\newcommand{\image}[2][0.5] {
		\begin{figure}[H]
			\centering
			\includegraphics[width=#1\paperwidth]{#2}
		\end{figure}
	}

	% optional parameters: height
	% required parameters: image1, label1, image2, label2, description
	\newcommand{\twoImages}[6][125px] {
		\begin{figure}[H]
			\centering
			\begin{subfigure}{.4\textwidth}
				\centering
				\includegraphics[height=#1]{#2}
				\caption{}
				\label{#3}
			\end{subfigure}
			\begin{subfigure}{.4\textwidth}
				\centering
				\includegraphics[height=#1]{#4}
				\caption{}
				\label{#5}
			\end{subfigure}
			\caption{#6}
		\end{figure}
	}

	\newcommand{\todo}[1]{
		\textcolor{red}{[TO-DO: #1]}
	}

	\newcommand{\eq}[1]
		{Eq.~#1}

	\newcommand{\eqs}[1]
		{Eqs.~#1}

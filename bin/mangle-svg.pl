#!/usr/bin/perl -w
#
while (<>) {
  next if /^\s+fill: none;$/;
  s/font-family: Liberation Sans;/font-family: Arial;/g;
  if (/<text x=.*Arial;.>(Node_\d+)<\/text/) {
    $id = $1;
#    s/$id//;
    s/<text x/<text id="$id" x/;
#  } elsif (/<text x=.*Arial;.>(.*)__(.*)__([YN])<\/text/) {
  } elsif (/<text x=.*Arial;.>(.*)__(.*)__(.*)<\/text/ ||
           /<text x=.*Glyphs.>(.*)__(.*)__(.*)<\/text/) {
    $id = $1;
    $treename = $2;
    $public = $3;
    s/<text x/<text id="$treename" x/;
#    if ($public eq "Y") {
    if ($public eq "Public") {
      s/>$id._$treename._$public/>$id/
    } else {
      s/>$id._$treename._$public/>/
    }
  }
  s/textLength='.*?px' lengthAdjust='spacingAndGlyphs'//;
  print;
}

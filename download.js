var data = source.data;
var headers = ["SMILES"];
for (header of Object.keys(data)) {
   if (header !== "IMGS" && header !=="SMILES") {headers.push(header)}
}

var filetext = headers.join("\t");
filetext = filetext.concat("\n");
for (i=0; i < data[headers[0]].length; i++) {
    var currRow = [];
    for (header of headers) {
        currRow.push(data[header][i].toString())}

    var joined = currRow.join("\t");
    joined = joined.concat("\n");
    filetext = filetext.concat(joined);
}

var filename = 'download.smi';
var blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename);
}

else {
    var link = document.createElement("a");
    link = document.createElement('a')
    link.href = URL.createObjectURL(blob);
    link.download = filename
    link.target = "_blank";
    link.style.visibility = 'visible';
    link.dispatchEvent(new MouseEvent('click'))
}

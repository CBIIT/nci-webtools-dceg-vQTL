import { Component, OnInit, ViewChild } from '@angular/core';
import { MatSort, MatPaginator, MatTableDataSource } from '@angular/material';
import { QTLsResultsService } from '../../../services/qtls-results.service';
import { environment } from '../../../../environments/environment' 

export interface Variant {
  gene_id: string;
  gene_symbol: string;
  variant_id: string;
  rsnum: string;
  chr: string;
  pos: string;
  ref: string;
  alt: string;
  tss_distance: number;
  pval_nominal: number;
  slope: string;
  slope_se: string;
  R2: string;
}

@Component({
  selector: 'app-qtls-locus-table',
  templateUrl: './qtls-locus-table.component.html',
  styleUrls: ['./qtls-locus-table.component.css']
})
export class QTLsLocusTableComponent implements OnInit {

  locusAlignmentData: Object;
  locusAlignmentDataQTopAnnot: Object;
  VARIANT_DATA: Variant[];
  selectedPop: string[];
  newSelectedPop: string;
  requestID: number;
  displayedColumns: string[] = ['gene_id', 'gene_symbol', 'variant_id', 'rsnum', 'chr', 'pos', 'ref', 'alt', 'tss_distance', 'pval_nominal', 'slope', 'slope_se', 'R2', 'LDpop', 'GWAS', 'genomAD'];
  dataSource = new MatTableDataSource<Variant>(this.VARIANT_DATA);

  @ViewChild(MatPaginator) paginator: MatPaginator;
  @ViewChild(MatSort) sort: MatSort;

  constructor(private data: QTLsResultsService) { }

  ngOnInit() {
    this.data.currentMainData.subscribe(mainData => {
      if (mainData) {
        this.locusAlignmentData = mainData["locus_alignment"]["data"][0]; // locus alignment data
        this.newSelectedPop = mainData["info"]["inputs"]["select_pop"][0]; // inputted populations
        this.locusAlignmentDataQTopAnnot = mainData["locus_alignment"]["top"][0][0]; // locus alignment Top Gene data
        this.requestID = mainData["info"]["inputs"]["request"][0]; // request id
      }
      this.selectedPop = this.newSelectedPop.split('+');; // recalculated new population selection
      this.expandPopulationGroup();
      if (this.locusAlignmentData) {
        this.VARIANT_DATA = this.populateVariantDataList(this.locusAlignmentData);
      }
      this.dataSource = new MatTableDataSource<Variant>(this.VARIANT_DATA);
      this.dataSource.paginator = this.paginator;
      this.dataSource.sort = this.sort;
    });
    // this.dataSource = new MatTableDataSource<Variant>(this.VARIANT_DATA);
    // this.dataSource.paginator = this.paginator;
    // this.dataSource.sort = this.sort;
  }

  populateVariantDataList(geneData) {
    var data = [];
    for (var i = 0; i < geneData.length; i++) {
      var variant = {};
      variant['gene_id'] = geneData[i]['gene_id'];
      variant['gene_symbol'] = geneData[i]['gene_symbol'];
      variant['variant_id'] = geneData[i]['variant_id'];
      variant['rsnum'] = geneData[i]['rsnum'];
      variant['chr'] = geneData[i]['chr'];
      variant['pos'] = geneData[i]['variant_id'].split(':')[1];
      variant['ref'] = geneData[i]['ref'];
      variant['alt'] = geneData[i]['alt'];
      variant['tss_distance'] = geneData[i]['tss_distance'];
      variant['pval_nominal'] = geneData[i]['pval_nominal'].toString();
      variant['slope'] = geneData[i]['slope'];
      variant['slope_se'] = geneData[i]['slope_se'];
      variant['R2'] = (geneData[i]['R2'] ? geneData[i]['R2'] : "NA").toString();
      variant['LDpop'] = "Go to";
      variant['GWAS'] = "Go to";
      variant['genomAD'] = "Go to";
      data.push(variant);
    }
    console.log(data);
    return data;
  }
  
  applyFilter(filterValue: string) {
    this.dataSource.filter = filterValue.trim().toLowerCase();
  }

  linkGTExGeneID(gene_id) {
    var url = "https://gtexportal.org/home/gene/" + gene_id
    var win = window.open(url, '_blank');
    win.focus();
  } 

  // linkGTExGeneSymbol(gene_symbol) {
  //   var url = "https://gtexportal.org/home/gene/" + gene_symbol
  //   var win = window.open(url, '_blank');
  //   win.focus();
  // } 

  linkGTExVariantID(chr, pos, ref, alt) {
    var url = "https://gtexportal.org/home/snp/" + chr + "_" + pos + "_" + ref + "_" + alt + "_b37"
    var win = window.open(url, '_blank');
    win.focus();
  } 

  // linkGTExRSNum(rsnum) {
  //   var url = "https://gtexportal.org/home/snp/" + rsnum
  //   var win = window.open(url, '_blank');
  //   win.focus();
  // } 

  linkLDpop(rsnum) {
    var QTopAnnotRef = this.locusAlignmentDataQTopAnnot["rsnum"];
    var selectedPopString = this.selectedPop.join('%2B');
    var url = "https://ldlink.nci.nih.gov/?tab=ldpop&var1=" + rsnum + "&var2=" + QTopAnnotRef + "&pop=" + selectedPopString + "&r2_d=r2"
    var win = window.open(url, '_blank');
    win.focus();
  } 

  linkGWAS(rsnum) {
    var url = "https://www.ebi.ac.uk/gwas/search?query=" + rsnum
    var win = window.open(url, '_blank');
    win.focus();
  }

  linkGenomADBrowser(chr, pos, ref, alt) {
    var url = "http://gnomad.broadinstitute.org/variant/" + chr + "-" + pos + "-" + ref + "-" + alt
    var win = window.open(url, '_blank');
    win.focus();
  }

  downloadTable() {
    var url = environment.endpoint + "tmp/" + this.requestID + ".variant_details.txt";
    var win = window.open(url, '_blank');
    win.focus();
  } 

  unique(value, index, self) {
    return self.indexOf(value) === index;
  }

  containsAll(subarr, arr) {
    for (var i = 0, len = subarr.length; i < len; i++) {
      if (!arr.includes(subarr[i])) {
        return false;
      }
    }
    return true;
  }

  remove(element, src) {
    var newArray = JSON.parse(JSON.stringify(src));
    // console.log(newArray);
    for (var i = 0; i < newArray.length; i++) {
      var idx = -1;
      if (newArray[i] == element) {
        idx = i;
      }
      if (idx != -1) {
        newArray.splice(idx, 1);
      }
    }
    return newArray;
  }

  expandPopulationGroup() {
    var african = ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"];
    var mixedAmerican = ["MXL", "PUR", "CLM", "PEL"];
    var eastAsian = ["CHB", "JPT", "CHS", "CDX", "KHV"];
    var european = ["CEU", "TSI", "FIN", "GBR", "IBS"];
    var southAsian = ["GIH", "PJL", "BEB", "STU", "ITU"];
    // AFR
    if (this.selectedPop.includes("AFR")) {
      this.selectedPop = this.remove("AFR", this.selectedPop);
      this.selectedPop = (this.selectedPop.concat(african)).filter(this.unique);
    }
    // AMR
    if (this.selectedPop.includes("AMR")) {
      this.selectedPop = this.remove("AMR", this.selectedPop);
      this.selectedPop = (this.selectedPop.concat(mixedAmerican)).filter(this.unique);
    }
    // EAS
    if (this.selectedPop.includes("EAS")) {
      this.selectedPop = this.remove("EAS", this.selectedPop);
      this.selectedPop = (this.selectedPop.concat(eastAsian)).filter(this.unique);
    }
    // EUR
    if (this.selectedPop.includes("EUR")) {
      this.selectedPop = this.remove("EUR", this.selectedPop);
      this.selectedPop = (this.selectedPop.concat(european)).filter(this.unique);
    }
    // SAS
    if (this.selectedPop.includes("SAS")) {
      this.selectedPop = this.remove("SAS", this.selectedPop);
      this.selectedPop = (this.selectedPop.concat(southAsian)).filter(this.unique);
    }
  }

}
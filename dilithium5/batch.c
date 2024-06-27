//
//  PQCgenKAT_sign.c
//
//  Created by Bassham, Lawrence E (Fed) on 8/29/17.
//  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
//
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "rng.h"  
#include "sign.h"
#include "speed_print.h"
#include "cpucycles.h"

#define	MAX_MARKER_LEN      50


#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR      -3
#define KAT_CRYPTO_FAILURE  -4

#define NTESTS 10000

uint64_t t[NTESTS];

void
fprintBstr(FILE *fp, char *s, unsigned char *a, unsigned long long l)
{
	unsigned long long  i;
 
	fprintf(fp, "%s", s);

	for ( i=0; i<l; i++ )
		fprintf(fp, "%02X", a[i]);

	if ( l == 0 )
		fprintf(fp, "00");

	fprintf(fp, "\n");
}
void hexStringToUint_8_a(const char* hexString, uint8_t* byteArray, size_t arrayLen) {
    
    for (size_t i = 0; i < arrayLen; i++) {
        sscanf(hexString + 2 * i, "%2hhX", &byteArray[i]);
    }
}

void messageGeneration(uint8_t** msgs, size_t mlens[]){
    const char* messages[] = {
        "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8",
        "225D5CE2CEAC61930A07503FB59F7C2F936A3E075481DA3CA299A80F8C5DF9223A073E7B90E02EBF98CA2227EBA38C1AB2568209E46DBA961869C6F83983B17DCD49",
        "2B8C4B0F29363EAEE469A7E33524538AA066AE98980EAA19D1F10593203DA2143B9E9E1973F7FF0E6C6AAA3C0B900E50D003412EFE96DEECE3046D8C46BC7709228789775ABDF56AED6416C90033780CB7A4984815DA1B14660DCF34AA34BF82CEBBCF",
        "2F7AF5B52A046471EFCD720C9384919BE05A61CDE8E8B01251C5AB885E820FD36ED9FF6FDF45783EC81A86728CBB74B426ADFF96123C08FAC2BC6C58A9C0DD71761292262C65F20DF47751F0831770A6BB7B3760BB7F5EFFFB6E11AC35F353A6F24400B80B287834E92C9CF0D3C949D6DCA31B0B94E0E3312E8BD02174B170C2CA9355FE",
        "1CDF0AE1124780A8FF00318F779A3B86B3504D059CA7AB3FE4D6EAE9FD46428D1DABB704C0735A8FE8708F409741017B723D9A304E54FDC5789A7B0748C2464B7308AC9665115644C569AE253D5205751342574C03346DDDC1950A6273546616B96D0C5ECE0A044AF0EDEFBE445F9AE37DA5AFB8D22A56D9FD1801425A0A276F48431D7AF039521E549551481391FE5F4EBFB7644D9F9782D83A95137E84EA3AEB3C2F8099",
        "DBE5B6C299B44F8D60FA972A336DF789EF4534EC9BA90DF92AD401D1907951EB6285EDA8F134277AB0A1145001C34E392187122506AA2DBB8617D7943A129EB5C07DF133D7CCDE94A7CB7F1795C62493ED375353D1F044257DA799F7D112C174FBC35687E2F87FEFBE2D83D29D7314B30A749FE41B1B81095638F112BC4563420AF235280E466FFBE7050C4937C60FC18D1A6025BCBD489F0C538E088E906ABE8597E2C8EBB64F01D225C847AAE4B77BAE6EBA9269962C4B94A9732CEAA2CB4093D442FFBCDD",
        "0073BEE97FC97C0FBC750D474AEB93189F061E1A5CF6600C04FB0464338EC7E85252F94FCBC7B2BD00E438480D9AF3ADD92A92E3E2E8ACB55077C3278FC7503988A76E9B6062996B20889AA55B343D5A003C8A8852D738F955799FA3426BE5CCD3AA6B6EDA04D4884941FFC0B69C5ACF12B347A74D0580CC3335BA816200F87674A4C1D98097C70F2F27C74E94A661850610ECF4847AB5B58344F958C5719E06BA396225BBE21ACB0FDC512B885D391E11B0C0ED5CE6B5DD8FAFF91F50025C69D43072F7706D80D9FD786E1104125D79A5F4B5FD838815D44FC8B1AB678078CC174DDE970D448B",
        "A1586245D81F96BD8EE81AA30F10C0ADB343D74CF72C4DFF71550C12873AF89FA1874D4731C996243C3749AF3F6188FFE9FA45430549045134EB29EF3CEC37E72904AA082B1C6161E6B52361E49AF4933A8D8C0734F21CAFD7467B0C02876F43211D6122E3E735FE36064DF7A0C91449237C2BC7C3A78AC7BB0F9567F2576F05802C872ADF183A87AA3B8217188F2F3535F877724F35B29E545DE4BCF258F13BBC7EDD8C6587F733C9691F74B4151CF8C060C3AE9E8D49FE7C77BF477DC9F23FD0F0B67320275529034B84F94176730923C03AA50F9584D9C2D60B8DCCF85A13F243F30A51ABEFBBF2CDA602BF3D75E849EB92422B808416C7E56B046CE38E4677AD24D23D7237A9",
        "9366ED7B3B623C411448B634446F1A3FAABDD163A6CC1E2BCAE4A98703CD8CEE441405892FBA051BE2A586A6950A5EF73A255E5F86B0D7212E0C51C3BC79BE4B88E76ED6F043FEF3204FAF044BFB1ED722D61EB5D0B74C66A257E8AC3A2206273C80D2EC2123A4DBB715D60118D99ED7322E38F1562F82379138DA3DDB8BAA7CE61AB729AFC3748C0134633CF45A9973C05C75D04E82F631845427626B5799DC07DDF830BA01E8BC6236BB6D03B37D949DBB29EEC7DFE60FBC17EA590956D251539792016E2A8B01E70476961BC9ADA43CDA682D0CAA4FCC58810BBA1A673EF8F6BC90BAEE701E8E4F7C04A346CA56C7B2862FF57756CE6CD1EE22D677BCDAA896EAE96F87870E032C18B6C6A0C1A191FAE2ED487CE55296CC4B6339EAC9E8A742BD0A44C3525CC750",
        "0998114C84F84080E7EEBB47D248980FAC9D28F1ABB6DBAB3DD59A5CFD2C7CFF7F308372874DD5447C7B02E30165501C0C673128E4C543A414222BDF47E7F4E8DCA757B0F4A3281C0D10C4F02AB52AAF5B9A715E012607BA310947A60A5F62D6B8CFA96386D27CFA709189202421C078934AA2D955468E550AD4D0D4ACDD98B168A9568E232192E92789830317FBC959087FFFE353B6C168F3EFBE7164444F1D6CBA5246E31658C65440A841DBA78257E78502843EC1A6E9710229C8EEB85D6CDDC7D543285624AA1F756A5DD4F1A5D4FA52DB8C5C34880ED448FBB6D254509FBEEA0FA022F276B6A66BEF7ABFEA6049FF74291BABE781F718683397077B29FA9E2B46BC6B09251E587CC5B182195DD4060CC4A319BFBE251A5B660A739DFE5D0E5B93F3CB7E440194F1C8BDA922CB1A3EE3D27EDFD61C1D31A7F4534E84889EC83B51F1641892766434",
        "4CCA95CB9F254C2EAA7DCFFEF662EE03320D5FC626A6484304BF62FC20F341FBE26E1537D7BD20E95440F7CC95EE84E1297C807A0BC9006DFCD5C22A5C1FC0865F5D70E5D63AD677FFFDEA52BF85D1A4F159F7ED16A745B4D971B620048B5F518EB2DC672CA35022578059E1ADAD7C07FE910A5D566B8321D9A12F34C250BE35CE964DDDEA23C90EA77C9C1BBE3532FEEFDA3637157786EC7D37775AE5CB0BB92EAB45A0FB1E833E8A6F3D06B85946E31A79B64A02B31FA640ED514A85882C89F693A06354DFDDB0B5E23E7792134C69C1D3908882DF3A7694A05B241B87FB2DBD1A4D9F26943B69F3CDF730301663089D1EBFC23299DA21300F735CEDF7B109F3E0BBE273776E6AAFA7054A6CD9682B967EB7903DE549E9558E62DCF3AC444DD7042FEA362EFB555BB97FB464AD7FAEABA3197C14A6740477DB50CE3FB8B762F48F880381D510FCC836E5880B48F08BD6333202E838AB73F2E106CFBFB218AAB802DA8A00F13F78FFB70C",
        "5C4B2E1A344DA1418B0F4BE3FD99505FC30F2A1E5B696E943BEE2451D7B268F722E04F8E00FDD9E1A470F8C977A6D45A5F621B8815E352FA14F64977D1FA08082A48AF495719EA6AC1C0B3D898603B4CF7EC88E68DD7190884382896D953D612CC21ABECFB01A04A1BB1BBE8986D34625756396CCD84BD1A6B5454DDA98824CD4844D98F356AB485EEB19F9196ABB1C3088C0C3C5846C88760B696D91A232D6F4CFFC85BFF33DE1A3433A27A209A461FCF37F2289F98BEA7CCF183DB1FC42A7EDF958E7913F8711DC375E43F09BE7C7A2C2B1318AE2A9CF5988FBC2CE0735A2CD9FB6C8496C34406C538C01BD494193240BFF947FED47B7CCE99A1747973F1FAA5223AC564BBA0CA8973D1310B5BFA1452CACE9110BC22A8D4080A8BAAA8ADFA3CFB6685679B648484E3A43F9B1B2531949BBB8FAE1846F6D45D9272FC2CAA2913B5D9F8D322E9B18A685122D74634C60730C101578BEF2480711FEFFE02123E76D6C846559E2EA99A98923EF095630102A5573EF027E0AB6E52555A9EDE0D15A73C8B2FEF87CA6FD9F903F0",
        "49755A7B1A7CDC5C9BDF5149968061D3C95EE67BFBAF02750C45094303A9D9CD23A08F19B9C768ADC63FFD1527186D09CA4E0356BB882E263BF015CBE3716C05B31A69DDDB790BA82C341AC9B6BE68A81B8BEF8D882304BAF0020D761A0DB04412033DC369961A5213B04E81736A580F1162780599CC029E262D67F31B2773AFB457A1ADAAA292163144F17DE384234F3303111FCD89BCB30333C6C6486F775ED099043C34E6C86450B650F1A02D03781B1D20691B767D166DADF1DCC4D8604D976EFDC9168373A7316DDA9B9FB02A4A321218D9F54E287B7167A08BC0153843BD6355AEA1310824DD5D5EC458BE694AF176119D9E588A29C650FF5500293659EA478B39A62149F819CDB7E7CB32E1D7B1284F159E2AB1B1EA41AF4D0AC94FF3111FC1CCD818F9B2CC7A259701405FDF6A51D2D3EF62789297BD16A659F14968EF902C4A23DA409BF13A4913467B5C991854B2CA6CC006D3F4197A6AA58BD5DD95C36928DA9583332C3FB134FA3890FE7E299F1C17205366C4F4230724C43E4803912E72B816658BBB1B63780865A1F66A2A49B96E93711B1BE97B827D12173402828B1A065B94310D5BD6098D",
        "439529DF1864297E33956AFEE00A60099B658A67830A6A6ABDDC329E87831D9F9B647917FEDF1AE182A40402143285516FCAB83F447354C72FAE81AC26E7005C2AA561763C152E66BD80F14565F47DEFA440DBB491E7994AB9FE35995D5FBB3800CA030B43DF611141637A5246AB9D9CAC02EFE14AF60736B6BDB2BABB97CF21E831E5D04D41C00F090B154977900EFADD3A9313389A3F84CB3AC38E8B57B70A43DD08A8243F8154013FD5CF29DE5A8DF0B197C12B17E0610FCFE3625CC94067E01E23D23A243AD1C1F805CC50E1447D1DF93C25B8D76396BB7199E64129522462C5FC8B30C132D4EE9E0BF6F52961FCE7ECF650647E7064AA5A6574649A323E144D7C5491DE4C0A1A76D08F93F87A2FC7F6955FEF86991E62E2CB42908E83B0C0A8BC180B7453CED293F1E20F300431EC1D395E8A537F0BC36A673D491F14381DEA90D8F176D06031B0A7AFB40EA8F76D37FA82E2572B9799A5FC7CF4C49BC20AD78EFA8CD989A84D72ED680AC3C0F64155C56ACBFD7C7D628B418A489F961357F77BD62204ADB079DD3106485A37FEE535C9CF82E832D8AADCBF686976B806B02AE733DB46DB0BF162E973931C3E338CC86DB38C66262D1B2EBC7691B8281E0B20BF36305FBA996D20ECFDC695",
        "8CB18850E27D8416B88A9A71F4A66BDF447814DB6C82098C371B53F61600EF5DFD88E4FB34200207C3F6F55166AF4878D38FCA7E2DC18FE662E3EA491B58A86246CAE16090FB7ADA53B9A67B3D0E3787D3323EA921274C60CFFB19A889BCF0300FE10E242AAE025F374DD83FBE9D007C8B9D9D75574C74146331DDEC6F0E49C10DBAF15654897E33E2B4780DBA484224AA6FAC79015D5792FAA2D532BB7D239B11D91420B98690B1FBDE9632223927E0804BFB284368A426C414C3DB8EA82F0D246413861475ED2DCA9E80FB4F3C34FEF7528069AE1975AFC52AC5AD2CDBCA1459E140F655556093210D7905A1A1E6CEEAEF0194A0B2EAB2C1EE853484E715D2A1DB551FDC620D5331164C74CA4848B61D408D2F2A943FA09EFEB63D524691C99DCC0B22CC61B98E6FB8039E5E0B2D7DE2CAAA900A44184BD56C9F02141A3AE8AFC661E3E898ECD3004FDB0704272BA780CD5DE35153B6FE223843024273642DCF8E4B58BE2AB1F61668680084AA0B75A32E766C8AE5EB30D4E02A12E6798DEA40F80D8DDFAD2041A52922701C689F46F49F84CFC05ECA6D7D4C356D50B6A0BA61966245D45134D6A1F5197540A1C39C36BB0B78831AF3F5156E669FD9213B64E0CF1C5A31E88AE79AD61757EC67B551B9F0A760F646BF81F6B92403A62840CC29FA4F3949B3A9F0A9A4286EE7808A",
        "9B64813C058F07A09A796FD764604EAF58CE144363702896DF0AB5FF26D5DE000D14BB8FD358FF5532D3B909AB62C18AC30F1900F84EBD3F4F18BD532D16C7B3470F0F8BDF72938C916DB18BCF1429DC1635B1C152C5F89A9EDB17116C11815A6C06273A889132923DA908FF39F4940A840D3CB575DC4D637AAFD37968EC61FC4EA04B4C320491A73ECFBDD8E10F1DFE902FCCEF93DD287ED872F67146BB8CA5A6ADCF0350E8BBA7F2F9762C4AA748FCE19748EB17334146C152FD63FAE3DFBB1A2C2B3C78960369551FDAC5D54643BEEAA59C1FEB0C21DBBB19977D848CD82A7AE0005F45956E0FE4700F14FBAA0C12FB8C65A6AEC95C5A5C8E79A6DA9C4E446872575C06AE49A31B82245E1757C7CE84D6D5DF3F642D3434B7E1A15A8B8A9DB460826B6CDCA69022DBF87595B582DDBB90A81E09A13C2AB1C125E4435FF30ABC9C56A00EDFA979F79D9C895E800D2DD6372FAE5FAACD83ADF8A6D55279D52DF547E9BAB39D99076AD7D297371344D35BD584E0FB5932F92FD5183B9250CD180FC645BEF6028C405B0EF35DAF783428173F1F2482AA1363640F66AF0FE8ECACC0DAB84ABD2A1FB53AF44445698CF1DDF4C2EA214DD339BE004E75BF76E95CA5C16981ABA5540689C1C1F1DAF4D0F89D62CCB3496340D61E7D5F5156FD3EDD02EDFEC8FCDD0B231697B0E66F4A3AAF46117532F5EE2CB4D2B3B82B0BEAE0A45A482CE9A976CC99AA82BEB0FE08CB68C4",
        "922320F7439E492F13C272A5738FF7122DD7A6B2832632E1F7A653FEF3B8639BCB9E84F482F22A948EA17DDE6958489593D2CB268BB52DF8ED612F2317BD6847D1622CF0532CB499ADC432233B93B6F7B1866B38975AC87859AC49F91E8D235846775F9E6E6D052339C741EF6178016EDB3D0B1E3F3536667B3EA2D489F88D254B8582421A31461374F465D7AD62E896BE0857134707A70477FABC09FE0A5CC3B3F32911F5FF3806B878205525AF69007F50535DF05C33AF3B0D00E297AC7EAA012E1D863DD5DD5FA47FB09467DBAD8BC42EDBAB42A9625BFDB9FE578343297506A3B71CDC8D5919955AF4605FCB0C7164D96A187AFF65D0F6210FEF2D11BA08D90C4458542BE72E084577BE9E451B8B6F4909884BCC5D25316ADCCD0925664D4D91C2E56433C1B68C632B0CA56D856DF1EDD5E113D1F026B30DAC4FD648A504F8F6809C701C97BCAC2B99286CEF5C1C923200B1BF6141EE1CFC51C5E14554BC02D7E058970254D2C02948360ABC4DFB439E66946A8AD615147BD8A6CB0886211E8B15DFF3C72B6F8908CE56BBC1B40E838103202E9F188D98E07555DB61778F895F76FBD838B6D14209D28EB393668924AC0E61072CBD9F93B864904FF4302DCEA131B2CA16BB04959ACEE096B1963CE07F59AB505FCC8D89FE08FC58751965F2F5CA753D76D58705652D3B1505E0F720EDE3142DE9776FFE4AA0C8A25E76C7A04843377C59F1002844E89189E22F621467B813A98BF07540A1649264F14A6844D65692617F7A4D93FA9A23829E256626",
        "576289D10AB03D5699EAC322D349F55C547101E4424BFA43BBBA3747B79F075AE1153A7A0AC8BB51D24FC46B7604E42EFE4343FA34AA4EB16D918F25E8A4D67C860CCA3F7480E1221ED3AE13A138F079FC252C6D7BEBC55CB81B86E74F339614BEBCF7E8F4440DF8678B01A4A41B3AFB1D112FE1C4C8D8C6BFE9D3EE2A335D477C60FBF43B2E5FFFE1546F5172EF51CFFB2A772E1575EAC79B24D49FD77F0BE351233E57EE6DCC7E2E29994873ABD434D34ACE83400C026E27E27888EA0BDD1BDE5A3E55AA8B5F2FEB57B8B0A96CD831906297C8169D04F15843A3249C50523CF56A4E19492EA16927DBA8759B88A99E0D20820E51FC9B6A6863115CF05C5BC3F4C869EB5A87124DF5DB102D737F3899CFAA5FEA4DD62DC4FEDB1AAFF67906ADAF8968020EFA5B10190F70E5F2C0F0457E4341BD449201D3A80AEB791254EC1C46DDCEBC3896C6DF702509BA62CD446D275806438EB4C03132B2E6BD01BD2F832D1D3C053C48C5A9DB1C4A22B130C4C9E96A2BF4C2A8F7DE0217A52D9AA5AEEE5E6A49708237EAB60B4019A51390C3EF10572A73D436875BB8D7D78543F96376E4BF3BCAABB92F89215E8D1093F3B287945708B5514BD7E62654D3BDF34B29009C64829A0CBF33C54D7AB0E81B81BDDA93028B341AB1DFF3D752DC4A1E5F9636A5C46E137EA35919D99E6571C5370C6E804BD2E2ABF566F035D65CF8F97E3E8F2ECAFA153BC6D8EC2831667A37FC96D1C2DA40BA84D0FB041DEF32AADAEF3F98CAFA957F6552F79D28A36B8BA20A9452671DE1BE8AF5D66714232507EDB9FF657F3D7E5FA7320FC0359A5F99280D446283BC",
        "021E9C06A2E4EF63D1A61958620C40016783879080D44311E04F2A446BCAEE5A486D17FF0F356BA70FF1C2B55BF957A59202903AE349878CB822E04275E0AFAABC0803BB6CDE3741E0BF9FCE0C5D5C814977474533DC63F9ED4F32AC3477A3EC9893EF55186728C85B03F4C2E61CA7733E1706766AEB8FEA80E233E8761B57FD5A3CEF700196674B34A3A55F68B3368B688FB1DDC976FF48BA6A98E2D66023F291A3C617A56CCBDB8732B8C34369ED11F4CCEA8FC8F673AD9FA0FD8990BEF70AF44C617FDFA096695D0C94EA8E17554F4461DC776DB2F416448B17680FE4D29B09E57603D8EBF55771AF84D8D4B9097302901C25CB6D73932E67C323D12C8ACB0E74CB89755F7EB3999D4EAB5E1B775E6B5C29D9733697030A26F3B93B3F286DB0F2DBDA71E1F103878063E77919D8892EB6A34F821B603ED4A898A9F30D00FEEF20985FEF1A7B7AF70DD29C269E88687F005D551EF05EB0603FD38745AED4F5BF4C2FC09F0604C98AE3A89E46BBFE907B87A1672DE547D651F035F392A8D4DB5E7260F43953028E312B95B9F25FFF2C0C579218390411D13D9A25F22DE4C7AA05FD11781DB08977160D48E02372C7D826F5CAC37D1A9B4230BE99A2D13CC2E9B2B17F0A1044EB9E0A2FBA376D35CDD2BC05F57DCE4BBC3BF07A09BCDE369929E6250EFDC61689466B040AEA376B09453A2C16813BBB685B54A225C49008BA6811E8BB5B3627F8C281244FDF5533216D126ED0E64FDABEC533424BFF77FE722CC438CA7587C19D965F0BF085D8692C27C5C84A9DEE53256D978948D89ABDF9842E0B765BE6A507D8630CBC5CA7FA0FBCA1CECC78D2E536AA7B2B902C4379777AC0920D69C57CC4E6032252BDE99E1A555E80D4",
        "7BEDAFEBABBBFB863CE496475F54E69A905AFA45899C3D7C16CFC73E31597D2404AE7014612E4CBFA238EFAF5B396B0B7435ADA5DE817E013188C280423C68924E1FA2A33CA56E6B85B7CCA7F00D3A6151F0629C1B92A13573320E0025863BBA7F3EEB987EE1B1A6230B10765DFC1FEEA498AE4B83521188E7503B506259103CEFB370E3651B06DD4F08013FF3AB9E2430626B0BD584232948462D85C0F82DA07B96FC65F62A43CD2F132D1A1D691C085980DAD8796CCE2FA0B268395EAC3DA2CC400F30F75BE87316216980CE213B48651DDB9E294F8CDB2CA05D3F2A507E4A03E2849AA8062918AFB5BCE9E4C3ABF2FFD4751DDDCF08AB09E36A29B830F3BAC6FEEBEA084575472E6F4B239AF89965A72954769A83E391DE467934237B07D8884A6B14CAD034FBF9BD7531D50D742E234E227E1A2DAF77A2FFACC579525134B15186D81AE6E5538871024BD2897475D6EE5B11BC51EDBB928D98475073785A75B331BF3D2297165AE6CF95C3A05F06DF747498462054F58A5AC736F96014B1A8CDB319D030D06DAD9CAB2B913F35FC392E1FC4B027CDBE775D64B04F1076A7C8F44C360745F98E87B84C18AB76F84F373F635AF4C8A87DF08DD4507899BAD892FF8CC1EE534D3277B5B82095628B84A7D5582149CF46C50AA963B56B4B91966B106B4B2EAA45D83A10993E8F933370AB29C6606B7CCFC41B21C6B99F2B9AC643E24300B350FA199EC10E64E4AF19181F78E8C43B2FA796241DC42CC8992BDFCDC39E7BC41BE68CDCE4FBC47C996DB42E8249EEDC146C216B514430C705FC939B9EEF677AD87F9CEE3398551FA0DAF774302324A410F4A4F4FC035CFBE960B38C390441E92D9E5624A8745976BC88FA538E398712361B77AD4CA5FF038D9F6CE157EB8A6137420D4E57018275DCEEBC4E480A5D"
    };

    for (int i = 0; i < N_MSG; i++) {
        mlens[i] = 33*(i+1);
        msgs[i] = malloc(mlens[i] * sizeof(uint8_t));
        if (msgs[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(1);
        }
        hexStringToUint_8_a(messages[i], msgs[i], mlens[i]);
    }
}



int main(){
    char                fn_req[32], fn_rsp[32];
    FILE                *fp_req, *fp_rsp;
    uint8_t             seed[48];
    uint8_t             msg[3300];
    uint8_t             entropy_input[48];
    uint8_t             *m;
    size_t              mlen, mlen1,smlen1;
    int                 count;
    int                 done;
    uint8_t             pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    int                 ret_val;

    poly p[154];
    // generate keypair
    if ( (ret_val = crypto_sign_keypair(pk, sk)) != 0) {
        printf("crypto_sign_keypair returned <%d>\n", ret_val);
        return KAT_CRYPTO_FAILURE;
     }

    uint8_t** sigs = malloc(sizeof(uint8_t*) * N_MSG);
    size_t* siglens = malloc(sizeof(size_t) * N_MSG);
    uint8_t** msgs= malloc(sizeof(uint8_t*) * N_MSG);
    uint8_t** msgs_verify= malloc(sizeof(uint8_t*) * N_MSG);
    size_t mlens[N_MSG];

    // generate 20 messages
    messageGeneration(msgs, mlens);

    for (int i = 0; i < N_MSG; ++i) {
        sigs[i] = malloc(sizeof(uint8_t) * (CRYPTO_BYTES + mlens[i]));
        size_t copy_position = CRYPTO_BYTES;
        size_t num_bytes_to_copy = mlens[i];
        memcpy(sigs[i] + copy_position, msgs[i], num_bytes_to_copy);
        siglens[i] = CRYPTO_BYTES + mlens[i];
    }

    // batch sign message for NTESTS observations
    for(int i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_signature_20_2(sigs, siglens, msgs, mlens, p, sk);
    }
    print_results("crypto_sign_signature_20_2:", t, NTESTS);

    // verify sigs original
    for(int i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        for (int i =0;i<N_MSG;i++){
            m = (uint8_t *)calloc(CRYPTO_BYTES + mlens[i], sizeof(uint8_t));
            if (ret_val = crypto_sign_open(m, &mlens[i], sigs[i], siglens[i], pk) != 0) {
                printf("crypto_sign_open returned <%d>\n", i);
                return KAT_CRYPTO_FAILURE;
            }
            free(m);
        }
    }
    print_results("crypto_sign_open:", t, NTESTS);

    // verify batch verification
    for(int i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        for (int i =0;i<N_MSG;i++){
            msgs_verify[i] = (uint8_t *)calloc(CRYPTO_BYTES + mlens[i], sizeof(uint8_t));
        }
        crypto_sign_open_batch(msgs_verify, mlens, sigs, siglens, pk, p);
        for (int i =0;i<N_MSG;i++){
            free(msgs_verify[i]);
        }
    }
    print_results("crypto_sign_open_batch:", t, NTESTS);

    // original sign message for NTESTS observations
    for(int i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        for (int k=0;k<N_MSG;k++){
            crypto_sign(sigs[k], &siglens[k], msgs[k], mlens[k], sk);
        }
    }
    print_results("crypto_sign:", t, NTESTS);

    for (int i = 0; i < N_MSG; i++) {
        free(msgs[i]);
		free(sigs[i]);
    }
    free(msgs_verify);
    free(msgs);
    free(siglens);
    free(sigs);
    return 0;
}